module RankCompV3

using Pkg
using Distributed
using SharedArrays
#using Base.Threads
using MultipleTesting
using DataFrames
using CategoricalArrays
using DelimitedFiles
using CSV
using Statistics
using StatsBase
using HypothesisTests
using Distributions
using Random
using LinearAlgebra
using FileIO


export 
	reoa, 
	McCullagh_test, 
    McNemar_exact_test, McNemar_test, McNemar_Bowker_test, 
	plot_result, plot_heatmap


include("$(joinpath(@__DIR__, "..", "code","plot.jl"))")

"""
pseudobulk_group(group_expr, n_pseudo, g_name)

Generate a matrix of pseudobulk profiles from `mat` which stores single-cell RNA profiles. Each column represents a cell's profile. Each pseudobulk profile is generated from `n_pseudo` (default: 10) single-cell profiles. `g_name` is a grouping label that will be added before the pseudobulk merged sample name.


# Examples
```jldoctest
julia> pseudobulk_group(DataFrame(rand(0:32, 10, 6),:auto), 3, String7("group1"))
10×3 DataFrame
 Row │ group1_x1  group1_x2  group1_x3
     │ Int64      Int64      Int64
─────┼─────────────────────────────────
   1 │        21         17         51
   2 │        35         38         61
   3 │        14         22         29
   4 │        23         11         36
   5 │        45         43         24
   6 │         5         27         49
   7 │        12         36         45
   8 │        24         18         17
   9 │        34         27         16
  10 │        30         42         29
```
"""
# pseudobulk
function pseudobulk_group(group_expr::DataFrame,
							n_pseudo::Int64,
							g_name::String7)
	r_group, c_group = size(group_expr )
	cp_group  = ceil(Int,  c_group/n_pseudo)
	cp_group > 1 || @info "WARN: too few profiles to generate $n_pseudo pseudo-bulk profiles for the 'group' group"
	it_group  = collect(Iterators.partition(sample(1:c_group, c_group, replace = false), cp_group)) # Random-shuffle, then partition
	group_expr  = reduce(hcat, [sum.(eachrow( group_expr[:, i])) for i in it_group ]) # Matrix r_group x n_pseudo
	group_expr  = DataFrame(group_expr,  :auto)
	rename!(group_expr,string.(g_name,"_",names(group_expr)))
	return group_expr
end

#   If two genes have the same expression value,
#   a random order is returned.
function is_greater(x::Number, y::Number)
	if abs(x - y) < 0.1
		return rand(Bool)
	else
		return x > y
	end
end

#   Return the lower limit for the count of the identical REOs
#   for being significanlly stable. 
function get_major_reo_lower_count(sample_size::Int,
				   pval_threshold::AbstractFloat=0.01)
	pval_min = pvalue(Binomial(sample_size), 0)
	if pval_min < pval_threshold
		return -findfirst(x-> pvalue(Binomial(sample_size), x) > pval_threshold, 
					0:floor(Int, sample_size/2)) + 2 + sample_size
    else
		@info "WARN: even if all samples have the identical REOs, it still cannot reach the required significance, $pval_threshold."
		@info "      The max significance level is $pval_min."
		return sample_size
	end
end

"""
    sum_reo(c_size, t_size, c_sign, t_sign, c_count, t_count)

Calculate the contribution of a REO of gene pair (a, b) to the 3x3 contigency table of gene a (not for gene b).

Return, for gene a 
```jldoctest
            Treatment
	       a < b   a ~ b  a > b
C	 a < b   n11    n12    n13
t	 a ~ b   n21    n22    n23
r	 a > b   n31    n32    n33
l
```
Only one cell is 1 for which both row and column conditions are met;

other cells are 0
"""
function sum_reo(
			      c_size::Int32,
			      t_size::Int32,
			      c_sign::Int32,
			      t_sign::Int32,
			     c_count::Int32,
			     t_count::Int32
			 )
	c_res = c_count >= c_sign ? 3 : ((c_size - c_count) >= c_sign ? 1 : 2)
	t_res = t_count >= t_sign ? 3 : ((t_size - t_count) >= t_sign ? 1 : 2)
	return c_res, t_res
end


function compute_pval( c_ctrl::Int32,
    c_treat::Int32,
     n_ctrl::Int32,
    n_treat::Int32,
      reo_t::AbstractMatrix{Int32}
   )
    n_ref, = size(reo_t)
	data_3x3 = [0 0 0;0 0 0;0 0 0]
	for i = 1:n_ref
		data_3x3[sum_reo(c_ctrl, c_treat, n_ctrl, n_treat, reo_t[i,:]...)...] += 1
	end
    pval, stat... = McCullagh_test(data_3x3)
    return pval, data_3x3, stat
end

"""
    compare_reos(ctrl, treat, n_ctrl, n_treat, ref_gene_vec)

Construct the contigency table for each gene and then test the association.

"""
function compare_reos(ctrl::AbstractMatrix,
		     treat::AbstractMatrix,
		    n_ctrl::Int32,
		   n_treat::Int32,
	  ref_gene_vec::BitVector            # 1-D Bool Array, 0 (false, other gene), 1 (true, ref gene)
    )
	r_ctrl,  c_ctrl  = size(ctrl)
	r_treat, c_treat = size(treat)
    n_ref  =  length(ref_gene_vec)
	r_ctrl == r_treat|| throw(DimensionMismatch("'Ctrl' and 'Treat' matrices do not have the equal number of rows."))
    r_ctrl == n_ref  || throw(DimensionMismatch("The number of rows of expr data != the length of ref gene list."))
	i_ref  = (1:r_ctrl)[ref_gene_vec]
	n_val  = 10+5  #5 for McCullagph_test, 3 for McNemar_exact_test, 2 for McNemar_Bowker_test and McNemar_test
	result = zeros(Float64, r_ctrl, n_val)
	@info "INFO: Number of execution threads: $(Threads.nthreads())"
	Threads.@threads for i in 1:r_ctrl
		#TODO: To reduce the time cost, consider subsample the ref gene list here. 
		max_ref = 3000
		local ref_tmp = i_ref[ i_ref .!= i ]
		if length(ref_tmp) > max_ref
			"If there are more than `max_ref` reference genes, sample `max_ref` of them randomly"
			ref_temp = sample(ref_tmp, max_ref, replace = false) 
		end
		local reo_mat = zeros(Int32, length(ref_tmp), 2)
		for jd in 1:length(ref_tmp)
			j = ref_tmp[jd]
			reo_mat[jd, 1] = Int32(sum(broadcast(is_greater,  ctrl[i,:],  ctrl[j,:]))) 
			reo_mat[jd, 2] = Int32(sum(broadcast(is_greater, treat[i,:], treat[j,:]))) 
		end
		pval, cont, stat = compute_pval(Int32(c_ctrl), Int32(c_treat), n_ctrl, n_treat, reo_mat) 
		result[i, :] = [pval 1 cont... stat...]
		if i%1000 == 0
			print(".")
		end
	end
	println(".")
	Δ = sort(result[:,12])
	# Filter out 10% most deviated  Δs before estimating the standard deviation
	se= std(Δ[round(Int, r_ctrl*0.05) : round(Int, r_ctrl*0.95) ])
	pval = pvalue.(Normal(0, se), result[:,12], tail=:both)
	padj = adjust(pval, BenjaminiHochberg())
	result[:, 1] = pval
	result[:, 2] = padj
	return result
end


"""
    McCullagh_test(mat::Matrix{Int})

Peform the McCullagh's test on a square contigency table `mat`.

McCullagh's test for a kxk contigency table;

see Biometrika, 1977, 64(3), 449-453.

Test case (Table 1 in p. 452) 

Example
```jldoctest
julia> mat = [43 8 3 0; 2 2 5 3; 1 0 7 2; 0 0 1 5]
julia> McCullagh_test(mat)
(0.005469174895116946, 1.4504988072997458, 1.502600073417028, 0.5221345956920705, 2.778017046308073)
```
Expected N matrix
```jldoctest
3×3 Matrix{Int64}:
 14   4  0
  4  12  3
  0   3  6
```
Expected R vector
```jldoctest
[11 11 5]'
```
Expected Δ1 = 1.45, Δ2 = 1.50, std. error = 0.53

"""
function McCullagh_test(mat::Matrix{Int})
	r, c = size(mat)
    r == c || throw(DimensionMismatch("input matrix 'mat' should be a square matrix."))
    #Symmetric matrix N, the first equation in p.450
    N = zeros(Int64, r -1, c -1)
	for i = 1:r-1
		for j = i:r-1
			N[i,j] = N[j, i] = sum(mat[1:i, j+1:end]) + sum(mat[j+1:end,1:i])
		end
	end
	n = diag(N)   # Extract the diagonal elements as a vector
	D = diagm(n)  # Construct a diagnonal matrix from a vector
	R = zeros(Int64, r-1)
	for i = 1:r-1
		R[i] = sum(mat[1:i,i+1:end])
	end
    #If N is a singular matrix, return 1.0 directly
	if abs(det(N)) <= eps()
        return 1.0, 0, 0, 0, 0
    else
        Ni = inv(N)
		ω2 = Ni*n   
		nu = 1.0/(n'*ω2)        # p.451
		ω1 = (D*ω2) .*nu
        Δ1 = ω1'*log.((R .+ 0.5)./(n - R .+ 0.5))
        Δ2 = log((0.5 .+ ω2' * R)/(0.5 .+ ω2' * (n .- R)))
        v1 = 4 * (1 + 0.25 * Δ1^2)*nu
        v2 = 4 * (1 + 0.25 * Δ2^2)*nu
        se = sqrt((v1 + v2) * 0.5)
        z1 = Δ1/se
		pval = pvalue(Normal(0, 1), z1, tail = :both)
		#println(ω1, ω2, Δ1, Δ2, se, z1, pval)
        return pval, Δ1, Δ2, se, z1
    end
end

"""
	McNemar_Bowker_test(mat::Matrix{Int}; continuity_correction = true)

Peform McNemar-Bowker's symmetry test for a square kxk contigency table;

calculate P value with Chi-squared distribution
"""
function McNemar_Bowker_test(mat::Matrix{Int};
                 continuity_correction = true)
    r, c = size(mat)
    r == c || throw(DimensionMismatch("input matrix 'mat' should be a square matrix."))
    para = (r*(r-1))>>1
    diff = mat .- mat'
    if continuity_correction & r == 2 & any(diff .!= 0)
        y = abs.(diff) .- 1
    else
        y = diff
    end
    x = mat .+ mat'
    if any([x[i,j] == 0 for i in 1:r-1 for j in i+1:r])
        return 1.0, 0
    else
        stat = sum(y[i,j]^2/x[i,j] for i in 1:r-1 for j in i+1:r)
        return pvalue(Chisq(para), stat, tail = :right), stat
    end
end

"""
    McNemar_test(n01::Int, n10::Int; continuity_correction = true)

Peform McNemar's test for a 2x2 contigency table;

calculate P value with Chi-squared distribution
"""
function McNemar_test(n01::Int, n10::Int;
              continuity_correction = true)
    dif = n01 - n10
    sum = n01 + n10
    if continuity_correction & dif != 0
        dif = abs(dif) - 1
    end
    if sum == 0
        pval = 1.0
    else
        stat = dif*dif/sum
        pval = pvalue(Chisq(1), stat, tail = :right)
    end
    return pval, stat
end

"""
    McNemar_exact_test(n01::Int, n10::Int)

Peform an exact McNemar's test for a 2x2 contigency table; 

calculate P value with Binomial distribution
"""
function McNemar_exact_test(n01::Int, n10::Int)
    n = max(n01, n10)
    N = n01 + n10
    dist = Binomial(N, 0.5)
    pval = pvalue(dist, n, tail = :both)
    return pval, N, n
end


"""
Perform differential expression analysis iteratively.

Calculate P values for each gene, then add the non-DEGs to the reference list
1) Calculate the pvals for each gene;

2) Update the ref_gene list by removing the DEGs from the ref_gene list and place non-DEGs as ref_gene;

3) If the number of ref_gene remains the same as the previous list, the iteration stops and returns. 
"""

# Group-specific DEGs
function identify_degs(
		    data::AbstractMatrix,
           group::AbstractVector, # Group information for each column in data
      gene_names::AbstractVector,
      # threshold::AbstractMatrix,   # Threshold for each group, 2xg
	    pval_reo::AbstractFloat,
        pval_deg::AbstractFloat,  # P-value threshold for DEGs
        padj_deg::AbstractFloat,  # FDR threshold for DEGs
	    ref_gene::BitVector, 
          n_iter::Int64,        # Threshold for iterations
	  	  n_conv::Int64         # Threshold for convergence
	)
	r, c = size(data)
	glen = length(group)
	glev = unique(group)   # Unique group levels
	gnum = length(glev)
	c == glen|| throw(DimensionMismatch("'data' and 'group' do not have compatiable sizes"))
	gnum > 1 || throw(DimensionMismatch("Only 1 level in 'group1, at least 2 levels!"))
	gind = [group .== i for i in glev]
	gsi1 = sum.(gind)
	gsi2 = c .- gsi1
	#TODO
	# threshold = ceil.(Int32, hcat(gsi1*0.9, gsi2*0.7)')
    threshold = get_major_reo_lower_count.(Matrix(hcat(gsi1,gsi2)') , pval_reo)
	R = falses(r*r, 9, gnum) # BitMatrix, 0 (false) or 1 (true), for each group
	@info "INFO: Number of execution threads: $(Threads.nthreads())"
	@info "INFO: Start constructing REO table"
	for i=1:r - 1
		i%100 == 0 && print(".")
		Threads.@threads for j=i+1:r
			# REO     1     2    3
			#  Ctrl: i<j  i~j  i>j
			#  Treat i<j  i~j  i>j
			local reo = broadcast(is_greater,  data[i,:],  data[j,:])
			local nre = [sum(reo[gind[k]]) for k in 1:gnum]
			local not = sum(nre) .- nre
			for k=1:gnum
				local ic = nre[k] >= threshold[1, k] ? 3 : ((gsi1[k] - nre[k]) >= threshold[1, k] ? 1 : 2)
				local it = not[k] >= threshold[2, k] ? 3 : ((gsi2[k] - not[k]) >= threshold[2, k] ? 1 : 2)
				# c\t i<j  i~j  i>j
				# i<j n11  n12  n13
				# i~j n21  n22  n23
				# i>j n31  n32  n33
				# R
				#  1    2    3   4    5    6   7    8    9
				# n11  n12  n13 n21  n22  n23 n31  n32  n33
				R[(i-1)*r+j, 3*(ic - 1) + it , k] = 1 # for gene i 
				R[(j-1)*r+i, 3*(3-ic) + (4-it), k] = 1# for gene j
				if gnum == 2
					break
				end
			end
		end
	end
	@info "\nINFO: Start identifying DEGs"
	res  = gene_names
	# outerloop for each group
	for k=1:gnum
		i_iter = 0
		result = zeros(Float64, r, 15) #10 +5 for McCullagph_test, +3 for McNemar_exact_test, +2 for McNemar_Bowker_test and McNemar_test
		ref_gene_vec = ref_gene
		while i_iter < n_iter
			# Construct contigency table For each gene
			Threads.@threads for i in 1:r
				local cont = sum(R[(((i-1)*r+1):(i*r))[ref_gene_vec], :, k], dims = 1)
				local pval, stat... = McCullagh_test(Matrix(reshape(cont, 3, 3)'))
				result[i, :] = [pval 1 cont... stat...]
				i%1000 == 0 && print(".")
			end
			print("\n")
			Δ = sort(result[:,12])
			# Filter out 10% most deviated  Δs before estimating the standard deviation
			se= std(Δ[round(Int, r*0.05) : round(Int, r*0.95) ])
			pval = pvalue.(Normal(0, se), result[:,12], tail=:both)
			padj = adjust(pval, BenjaminiHochberg())
			# padj = adjust(result[:,1], BenjaminiHochberg())
			result[:, 1] = pval
			result[:, 2] = padj
			inds = .!((pval .<= pval_deg) .&& (padj .<= padj_deg))
			@info "INFO: iteration $i_iter,  # DEGs $(length(inds) - sum(inds)), # non-DEGs $(sum(inds))"
			if (abs(sum(ref_gene_vec) - sum(inds)) < n_conv)
				@info "INFO: Convergence threshold is reached"
				break
			end
			i_iter += 1
			ref_gene_vec = inds
		end
		gene_up_down = copy(gene_names)
		gene_up_down[result[:,15] .> 0] .= "up"
		gene_up_down[result[:,15] .< 0] .= "down"
		gene_up_down[result[:,15] .== 0] .= "no change"
		res = hcat(res, result, gene_up_down)
		if gnum==2
			@info "INFO: The results of differentially expressed genes in the iteration process of $(glev[1]) vs $(glev[2]) were output."
			break
		end
		@info "INFO: The results of differentially expressed genes in the iteration process of $(glev[k]) vs other were output."
	end
	return res
end

"""
    reoa(fn_expr::AbstractString = "fn_expr.txt",
           fn_meta::AbstractString = "fn_meta.txt")

RankCompV3 is a differential expression analysis algorithm based on relative expression ordering (REO) of gene pairs. It can be applied to bulk or single-cell RNA-sequencing (scRNA-seq) data, microarray gene expression profiles and proteomics profiles, etc. When applied in scRNA-seq data, it can run in single-cell mode or pseudo-bulk mode. The pseudo-bulk mode is expected to improve the accuracy while decreasing runntime and memory cost. 

# Examples
## Quick Start

Run a test job with the input files distributed with the package.

```jldoctest
julia> using RankCompV3
# Use the default values for the following other parameters. If you need to modify the parameters, add them directly.
julia> result = reoa(use_testdata="yes")
```

The analysis results and a few plots will be generated and saved in the current work directory. They are also returned by the `reoa` function and can be captured by assign the returned values to a variable,  e.g., `result` in the above example.  

The first return value is a DataFrame, where rows are genes and columns are statistical values for each gene. All the genes passing the basic preprocessing step are retained. 


```jldoctest
julia> result
(19999×16 DataFrame
   Row │ Name     pval         padj        n11      n21      n31      n12      n22      n32      n13      n23      n33      Δ1           Δ2          se         z1
       │ String   Float64      Float64     Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64      Float64     Float64    Float64
───────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │ DE1      0.23566      0.716036     1532.0     75.0      0.0   1220.0  11010.0    602.0     16.0   2933.0   1674.0   1.91208      1.81954    0.0393608    48.5783
     2 │ DE2      0.280118     0.761567     2277.0    288.0      0.0    965.0  11514.0    372.0     20.0   2615.0   1011.0   1.74141      1.70287    0.0405313    42.9647
     3 │ DE3      0.0578546    0.359121     1576.0     37.0      0.0   1376.0   8716.0    293.0    113.0   5211.0   1740.0   3.05828      3.02011    0.0437133    69.9623
   ⋮   │    ⋮          ⋮           ⋮          ⋮        ⋮        ⋮        ⋮        ⋮        ⋮        ⋮        ⋮        ⋮          ⋮           ⋮           ⋮           ⋮
 19998 │ EE19999  0.344847     0.823375     1475.0    307.0      0.0   1756.0  15356.0    128.0      0.0     38.0      2.0   1.52307      1.41599    0.0525798    28.9668
 19999 │ EE20000  0.694484     0.980397     1571.0    315.0      0.0    979.0  15555.0    362.0      1.0    225.0     54.0   0.63329      0.577392   0.0481845    13.143
                                                                                                                                             19965 rows omitted, 19999×16
```

## Run your own DEG analysis

You need to prepare two input files before the analysis: metadata file and expression matrix. Both of them should be saved in the `TSV` or 'CSV` format and they should be compatible with each other.   

- **metadata file (required).**

 The metadata file contains at least two columns. The first column is the sample names, and the second column is the grouping information. Only two groups are supported at present, therefore, do not include more than two groups. 

 Column names for a metadata should be `Name` and `Group`. 

 See an example metadata file, [fn_metadata.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_metadata.txt)


- **expression matrix file (required).**

 The first column is the gene name and the column header should be `Name` and the rest columns are profiles for each cell or each sample. Each column header should be the sample name which appears in the metadata file.

 See an example expression matrix file, [fn_expr.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_expr.txt)

Once the files are ready, you can carry out the DEG analysis with the default settings as follows. 

```jldoctest
julia> using RankCompV3
# Use the default values for the following other parameters. If you want to modify the parameters, add them directly.
julia> reoa("/public/yanj/data/fn_expr.txt",
		"/public/yanj/data/fn_metadata.txt")
```

Other parameters can be set by passing the value to the corresponding keyword. 

```jldoctest
julia> reoa("/public/yanj/data/fn_expr.txt",
    	"/public/yanj/data/fn_metadata.txt",
    	expr_threshold = 0,
    	min_profiles = 0,
    	min_features = 0,
    	pval_reo = 0.01,
     	pval_deg = 1.0,
     	padj_deg = 0.05,
    	use_pseudobulk = 0,
    	use_hk_genes = "yes"
    	hk_file = "HK_genes_info.tsv",
    	gene_name_type = "ENSEMBL",
    	ref_gene_max = 3000,
    	ref_gene_min = 100
    	n_iter = 128,
    	n_conv = 5,
    	work_dir = "./",
    	use_testdata = "no")
```

"""
function reoa(
	 	   fn_expr::AbstractString = "fn_expr.txt",
           fn_meta::AbstractString = "fn_meta.txt";
    expr_threshold::Number = 0,
	  min_profiles::Int = 0, # Include features (genes) detected in at least this many cells
      min_features::Int = 0, # Include profiles (cells) where at least this many features are detected
          pval_reo::AbstractFloat = 0.01,
          pval_deg::AbstractFloat = 1.0,
          padj_deg::AbstractFloat = 0.05,
		  n_pseudo::Int = 0, # 0 for not using pseudobulk mode, Other values indicate the number of samples that are combined after each pseudobulk.
      use_hk_genes::AbstractString = "yes",
           hk_file::AbstractString = "$(joinpath(@__DIR__, "..", "hk_gene_file", "HK_genes_info.tsv"))",
    gene_name_type::AbstractString = "ENSEMBL", # Available choices: ENSEMBL, Symbol, ENTREZID ...
      ref_gene_max::Int = 3000,    # If the number of available features is higher than this, take a random sample of this size
      ref_gene_min::Int = 100,     # If the number is lower than this, ignore the house-keeping genes
	        n_iter::Int = 128,     # Max iterations 
	        n_conv::Int = 5,       # Convergence condition: max. difference in the number of DEGs between two consective iterations
          work_dir::AbstractString = "./",
      use_testdata::AbstractString = "no"
    )

    cd(work_dir)
	# if use test data, we will ignore the input for 'fn_expr' and 'fn_meta'.
    if use_testdata == "yes"
        fn_expr = "$(joinpath(@__DIR__, "..", "test", "fn_expr.txt"))"
        fn_meta = "$(joinpath(@__DIR__, "..", "test", "fn_meta.txt"))"
    end

    # Import datasets
	isfile(fn_expr) && isfile(fn_meta)             || throw(ArgumentError("$fn_expr, or $fn_meta, does not exist or is not a regular file."))
	filesize(fn_expr) > 0 && filesize(fn_meta) > 0 || throw(ArgumentError("$fn_expr, or $fn_meta, has size 0."))
    fn_stem, = splitext(basename(fn_expr))   #filename stem
	# Read in expression matrix
    @time expr = occursin(".rds",fn_expr) ? DataFrame(load(fn_expr),:auto) : (occursin(".RData",fn_expr) ? DataFrame(get(load(fn_expr),fn_stem,1),:auto) : CSV.read(fn_expr, DataFrame))
    @time meta = occursin(".rds",fn_meta) ? DataFrame(load(fn_meta),:auto) : (occursin(".RData",fn_meta) ? DataFrame(get(load(fn_meta),splitext(basename(fn_meta))[1],1),:auto) : CSV.read(fn_meta, DataFrame))
	mr,mc= size(meta)
    # r,c  = size(expr)
	# Check the compatability between meta data and expression matrix
	mc >= 2 || throw(ArgumentError("$fn_meta the file for meta data, has only 0 or 1 column."))
	if !(["Name" "Group"] ⊆ names(meta)) && (meta[:,1] ⊆ names(expr))
		# assume the first two column as 'Name' and 'Group'
		rename!(meta, [1 => :Name, 2 => :Group]) 
	end
	# exit early if meta data does not fit the expression profile
	if !(["Name" "Group"] ⊆ names(meta)) || !(meta.Name ⊆ names(expr))
		throw(ArgumentError("Meta data file, $fn_meta, does not fit with the expression file, $fn_expr. Some sample names in the meta are not found in the column names of the expression matrix"))
	end
	if size(unique(names(expr)))[1] .!= size(names(expr))[1]
		throw(ArgumentError("Duplicate column names exist in the representation matrix."))
	end
	meta.Group = categorical(meta.Group)
	g_name = unique(meta.Group)
	mg         = length(g_name)
	if mg < 2
		throw(ArgumentError("Meta data file, $fn_meta has only 0 or 1 group. It must consist of two 'Group' levels"))
	end
	@info "INFO: According to the meta information, there are $mg groups of data and each group will be analyzed with the rest of the sample."
	# if expr has no column named 'Name', we check if the first column is not a data column (the column name appears in meta.Name)
	# then we consider the first column as the gene names 
	if !(["Name"] ⊆ names(expr)) && !(names(expr, 1) ⊆ meta.Name)
		"Expression matrix has no 'Name' column and its first column is not a sample profile; designate it as gene name column."
		rename!(expr, [1 => :Name])
	end
	# Remove the rows with missing data
	# Columns with Number datatype are assumed to be the expression profiles
	data_names = names(dropmissing!(expr), Number)
	if !(meta.Name ⊆ data_names)
		throw(ArgumentError("$fn_expr expression matrix contains non-numeric (Number) profiles."))
	end
    gene_names = convert(Vector{String},expr.Name)
	# Generate pseudo-bulk profiles
	if n_pseudo > 0
		df_expr = reduce(hcat, [pseudobulk_group(expr[:, meta.Name[meta.Group .== g_name[i]]], n_pseudo, g_name[i]) for i in 1:mg ])
		meta_group = DataFrame()
		insertcols!(meta_group,   1, :Name => names(df_expr))
		insertcols!(meta_group,   2, :Group => reduce(hcat,split.(names(df_expr),"_"))[1,:])
	else
		df_expr = expr[:,2:end]
		meta_group = meta
	end
	# Filter low-expressed cells
	s_inds = (sum.(eachcol(df_expr .> 0)) .> min_profiles)
	s_del = sum(.!s_inds)
	if s_del > 0
		for i in 1:s_del
			meta_group = meta_group[.!(meta_group[:,1] .== names(df_expr)[.!s_inds][i]),:]
		end
	end
	df_expr    = df_expr[:, s_inds]
	inds       = (sum.(eachrow(df_expr .> 0)) .> min_features)
	gene_names = gene_names[inds]
	df_expr    =  df_expr[inds, :]
	@info "INFO: size after filtering lowly expressed genes and profiles and pseudo-bulk sampling, $(size(df_expr))"
	# Process the house-keeping genes
    if gene_name_type == "ENTREZID"
        gene_names = convert(Vector{String}, gene_names)
    end
	# If not using house-keeping genes, randomly pick 'ref_gene_max' as reference genes
	ref_gene = sample(gene_names, min(length(gene_names), ref_gene_max), replace = false)
    if use_hk_genes == "yes"
        # Read in the list of house-keeping genes
        isfile(hk_file)       || throw(ArgumentError("$hk_file does not exist or is not a regular file."))
        filesize(hk_file) > 0 || throw(ArgumentError("$hk_file for house-keeping genes has size 0."))
		# Read in as String
        ref_gene = CSV.read(hk_file, DataFrame, types = String)
		if gene_name_type ∈ names(ref_gene)
			ref_gene = ref_gene[:, gene_name_type]
			ref_gene = intersect(ref_gene, gene_names)
			if length(ref_gene) < ref_gene_min
				@info "WARN: only $(length(ref_gene)) house-keeping genes are available, we just ignore this."
				ref_gene = sample(gene_names, min(length(gene_names), ref_gene_max), replace = false)
			end
		end
	end
	ref_gene_vec = convert(BitVector, [i ∈ ref_gene for i in gene_names])
	@time result = identify_degs(Matrix(df_expr),
								meta_group.Group, # Group information for each column in data
								gene_names,
								# threshold::AbstractMatrix,   # Threshold for each group, 2xg
								pval_reo,
								pval_deg,  # P-value threshold for DEGs
								padj_deg,  # FDR threshold for DEGs
								ref_gene_vec, 
								n_iter,        # Threshold for iterations
								n_conv         # Threshold for convergence
						)
	# Beautify output
	# McCullagh_test result header
	header = [:pval, :padj, :n11, :n12, :n13, :n21, :n22, :n23, :n31, :n32, :n33, :Δ1, :Δ2, :se, :z1, :up_down]
	for i in 1:mg
		g_sit = (meta_group.Group .== g_name[i])
		outcome_result = DataFrame(result[:,(2 + 16*(i - 1)):(1 + 16*i)],header)
		insertcols!(outcome_result,   1, :genename => result[:,1])
		(mg == 2) ? fg_name = join([g_name[1],g_name[2]],"_") : fg_name = g_name[i]
		CSV.write(join([fn_stem,fg_name,"result.tsv"], "_"), outcome_result,  delim = '\t')
		plot_result(outcome_result, join([fn_stem,fg_name], "_"))
		plot_heatmap(df_expr[:,g_sit], df_expr[:,.!g_sit], join([fn_stem,fg_name], "_"), log1p = true)
		n_degs = ((outcome_result.padj .<= padj_deg) .&& (outcome_result.pval .<= pval_deg))
		(sum(n_degs) != 0) ? plot_heatmap(df_expr[n_degs,g_sit], df_expr[n_degs,.!g_sit], join([fn_stem, fg_name, "degs"], "_"), log1p = true) : "no degs"
		(mg == 2) ? break : continue
	end
	insertcols!(df_expr,   1, :genename => result[:,1])
	CSV.write(join([fn_stem, "df_expr.tsv"], "_"), df_expr, delim = '\t')
	CSV.write(join([fn_stem, "df_meta.tsv"], "_"), meta_group, delim = '\t')
	@info "INFO: The expression profile and metadata file after pseudobulk are $(join([fn_stem, "df_expr.tsv"], "_")) and $(join([fn_stem, "df_meta.tsv"], "_"))"
	gene_up_down = DataFrame(hcat(gene_names, reduce(hcat, [result[:,i*16 + 1] for i in 1:fld(size(result)[2],16)])),:auto)
	CSV.write(join([fn_stem, "gene_up_down.tsv"], "_"), rename!(gene_up_down,((mg == 2) ? ["gene_name", string(g_name[1], "_vs_", g_name[2])] : ["gene_name"; string.(g_name, "_vs_other")])), delim = '\t')
    return gene_up_down
end

end # module
