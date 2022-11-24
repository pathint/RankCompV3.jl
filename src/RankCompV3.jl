module RankCompV3

using Pkg
using Distributed
using SharedArrays
using Base.Threads
using MultipleTesting
using DataFrames
using DelimitedFiles
using CSV
using Statistics
using RCall
using HypothesisTests
using Distributions
using Random
using LinearAlgebra
using ArgParse


export reoa

# Determine whether the required R package exists
function check_R_packages()
    R"""
    if(require("ggplot2")){
      print("ggplot2 package loaded successfully.")
    } else {
      print("No ggplot2 package exists, trying to install.")
      install.packages("ggplot2")
      if(require("ggplot2")){
        print("ggplot2 installed successfully and loaded package.")
      } else {
        stop("Installation failed.")
      }
    }
    if(require("pheatmap")){
      print("pheatmap package loaded successfully.")
    } else {
      print("No heatmappackage exists, trying to install.")
      install.packages("pheatmap")
      if(require("pheatmap")){
        print("pheatmap installed successfully and loaded package.")
      } else {
        stop("Installation failed.")
      }
    }
    """
    return 
end


#   If two genes have the same expression value,
#   a random order is returned.
function is_greater(x::Number, y::Number)
	if abs(x - y) < 0.5
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
	#floor向下取整
    else
		println("WARN: even if all samples have the identical REOs, it still cannot reach the required significance, $pval_threshold.")
		println("      The max. significance is $pval_min.")
		return sample_size
	end
end

function sum_reo(
			      c_size::Int32,
			      t_size::Int32,
			      c_sign::Int32,
			      t_sign::Int32,
			     c_count::Int32,
			     t_count::Int32
			 )
	"""
	The triple table
	Return 
	1 = (0, 0，1), a < b in ctrl and a > b in treat
	2 = (0, 1，0), a ~ b in ctrl and a ~ b in treat
    3 = (1, 0，0), a > b in ctrl and a < b in treat
	and so on
	"""
    if c_count >= c_sign
        #A>B,c
		if (t_size - t_count) >= t_sign
            # A>B,c;A<B,t
			return ([0, 0, 0], [0, 0, 0], [1, 0, 0])
        else
            if t_count >= t_sign
                #A>B,c;A>B,t
                return ([0, 0, 0], [0, 0, 0], [0, 0, 1])
            else
                #A>B,c;A~B,t
                return ([0, 0, 0], [0, 0, 0], [0, 1, 0])
            end
		end
	else
        if (c_size - c_count) >= c_sign
            #A<B,c
            if t_count >= t_sign
                #A<B,c;A>B,t
                return ([0, 0, 1], [0, 0, 0], [0, 0, 0])
            else
                if (t_size - t_count) >= t_sign
                    #A<B,c;A<B,t
                    return ([1, 0, 0], [0, 0, 0], [0, 0, 0])
                else
                    #A<B,c;A~B,t
                    return ([0, 1, 0], [0, 0, 0], [0, 0, 0])
                end
            end
        else
            if t_count >= t_sign
                #A~B,c;A>B,t
                return ([0, 0, 0], [0, 0, 1], [0, 0, 0])
            else
                if (t_size - t_count) < t_sign
                    #A~B,c;A<B,t
                    return ([0, 0, 0], [1, 0, 0], [0, 0, 0])
                else
                    #A~B,c;A~B,t
                    return ([0, 0, 0], [0, 1, 0], [0, 0, 0])
                end
            end
        end
	end
end

function compute_pval( c_ctrl::Int32,
    c_treat::Int32,
     n_ctrl::Int32,
    n_treat::Int32,
      reo_t::AbstractMatrix{Int32}
   )
    n_ref, = size(reo_t)
    obs_all = reduce(.+, map(x -> sum_reo(c_ctrl, c_treat, n_ctrl, n_treat, x...), eachrow(reo_t)))
    data_3x3 = [obs_all[1]';obs_all[2]';obs_all[3]']
    delta = McCULLACH_test(data_3x3)
    f = open( "delta_nd_z.tsv", "a" )
    m = vcat([delta[1]],[delta[2]])
    writedlm( f,[m] )
	close(f)
    delta = delta[1]
    return delta
end

function compare_reos(ctrl::AbstractMatrix,
		     treat::AbstractMatrix,
		    n_ctrl::Int32,
		   n_treat::Int32,
		  ref_gene::Vector{Bool},            # 1-D Bool Array, 0 (false, other gene), 1 (true, ref gene)
      ref_gene_new::BitVector
    )
	"""
	Compute p-value for each gene.
	"""
	 r_ctrl,  c_ctrl = size(ctrl)
	r_treat, c_treat = size(treat)
	r_ctrl == r_treat || throw(DimensionMismatch("the number of rows of 'ctrl' not equal to the number of rows of 'treat'."))
    r_ctrl == length(ref_gene)|| throw(DimensionMismatch("the number of rows of 'ctrl' not equal to the length of 'ref_gene'."))
    n_ref = sum(ref_gene)
    n_rep = c_ctrl + c_treat
	println("INFO: Number of threads, $(Threads.nthreads())")
	i_deg = (1:r_ctrl)[.!ref_gene]
	i_ref = (1:r_ctrl)[  (ref_gene_new .| ref_gene)]
    data = ["delta" "nz_d"]
    writedlm("delta_nd_z.tsv", data, '\t')
    pval_t1 = pmap(i -> compute_pval(Int32(c_ctrl), 
                       Int32(c_treat), 
                       n_ctrl, n_treat, 
                       reduce(vcat, [[Int32(sum(broadcast(is_greater,  ctrl[i,:],
                                        ctrl[j,:])))	Int32(sum(broadcast(is_greater, treat[i,:],
                                        treat[j,:])))] for j in i_ref[(i_ref .!= i)]])),
                i_deg)
    data = CSV.read("delta_nd_z.tsv", DataFrame)
    pval_t = data[:,1]
    sd_delta = data[:,2]
    pval_t_retain = (pval_t .!= 111)
    pval_t_del = pval_t[pval_t_retain]
    sd_delta = sd_delta[pval_t_retain]
    dist = Normal(0,1)
    pval_t_var = pvalue.(dist, sd_delta, tail = :right)
    padj = adjust(pval_t_var, BenjaminiHochberg())
    pval_t_del = convert(Vector{Float64}, pval_t_del)
	return pval_t_del,sd_delta,pval_t_var,padj,pval_t_retain
end


function McCULLACH_test(data::Matrix)
	l_i = size(data)[1]
    l_j = size(data)[2]
    #N,n,r
    N = zeros(l_i-1,l_j-1)
    D = zeros(l_i-1,l_j-1)
    n = zeros(l_i-1,1)
    r = zeros(l_i-1,1)
    m = 0
    for i in 1:l_i-1
        N[i,i] = n[i,1] = D[i,i] = sum(data[1:i,i+1:l_i]) + sum(data[i+1:l_i,1:i])
        r[i,1] = sum(data[1:i,i+1:l_i])
        for j in i+1:l_j-1
            if i != j
                N[i,j] = N[j,i] = N[i,i] - data[j,i] - data[i,j] - m
                m = data[j,i] + data[i,j]
            end
        end
    end
    #If the determinant is a singular matrix, mark 111.0.
    if det(N) == 0
        return 111.0,111.0
    else
        N_inv = inv(N)
        #w1,w2,delta1,delta2
        w1 = (D*N_inv*n)/(n'*N_inv*n)
        w2 = N_inv*n
        delta1 = w1'*log.((r .+ 1/2)./(n - r .+ 1/2))
        delta2 = log((1/2 .+ w2' * r)/(1/2 .+ w2' * (n - r)))
        var_delta1 = (4 * (1 .+ (1/4) * delta1^2)/(n'*N_inv*n))
        var_delta2 = (4 * (1 .+ (1/4) * delta2^2)/(n'*N_inv*n))
        sd = sqrt((var_delta1 + var_delta2) /2)
        nd_z1 = delta1/sd
        return delta1[1],nd_z1[1]
    end
end

function rand_sit(ref_gene_num::Int64,ref_sum_yuan::Int64)
    rand_num = unique([rand(1:ref_sum_yuan) for r in 1:ref_gene_num*100])
    if length(rand_num) >= ref_gene_num
        return rand_num[1:ref_gene_num]
    else
        return rand_sit(ref_gene_num)
    end
end

function hk_non_hk_graph(data::DataFrame,sample_name::String,fn_stem::String)
    @rput data
    @rput sample_name
    @rput fn_stem
    R"""
    library(ggplot2)
    getmode <- function(v) {
       uniqv <- unique(v)
       uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    sample_data=data[,which(colnames(data)==sample_name)]
    max_num=getmode(sample_data)
    hk_graph = ggplot(data, aes_string(x=sample_name, fill="gene_types"))+
    geom_density(alpha=.3)+
    theme_classic()+
    scale_x_continuous(limits = c(max(max_num-20,0),min(max_num+20,max(sample_data))))
    ggsave(hk_graph,filename = paste(fn_stem,"_hk_nonhk_gene_",sample_name,".pdf", sep = ""),width = 6,height = 5)
    """
end

function delta_p_padj_graph(data::DataFrame,fn_stem::String)
    @rput data
    @rput fn_stem
    R"""
    library(ggplot2)
    delta_graph=ggplot(data, aes(x=delta)) + 
      geom_histogram(aes(y=..density..),
                     binwidth=1,
                     colour="black", fill="white") +
      geom_density(alpha=.2, fill="#4b5cc466")+
      geom_function(fun = dnorm, args = list(mean = mean(data$delta), sd = sqrt(var(data$delta))),colour = "red")+
      theme_classic()
    ggsave(delta_graph,filename = paste(fn_stem, "_delta_graph.pdf", sep = ""),width = 6,height = 5)
    graph_p_padj=function(data,name){
      num = data[,which(names(data)==name)]
      bin_size = max(num)/100
      ts = matrix(seq(min(num),max(num)+max(num)/100,bin_size),ncol = 1,byrow = TRUE)
      ts = cbind(ts,ts)
      colnames(ts)=c(name,"num")
      nu=0
      for (i in 1:length(ts[,1])) {
        s=data[ which(num <= ts[i,1]), ]
        if (i==1){
          ts[i,2] = length(s[,1])
          nu=nu+ts[i,2]
          next
        }
        else{
          ts[i,2] = length(s[,1]) - nu
          nu=nu+ts[i,2]
        }
      }
      ts = data.frame(ts)
      ts[,2]=ts[,2]/sum(ts[,2])
      name_graph = ggplot(ts,aes_string(x = name, y = "num"))+
        geom_point()+
        theme_classic()
        ggsave(name_graph,filename = paste(fn_stem, "_",name,"_graph.pdf", sep = ""),width = 6,height = 5)
    }
    graph_p_padj(data,'Pval')
    graph_p_padj(data,'Padj')
    graph_p_padj(data,'sd_delta')
    """
end

function deg_exp_graph(data::DataFrame,fn_stem::String,c_ctrl::Int64,c_treat::Int64)
    @rput data
    @rput fn_stem
    @rput c_ctrl
    @rput c_treat
    R"""
    rownames(data)=data[,1]
    data=data[,-1]
    data=log2(data)
    data[data==Inf]<-0
    data[data==-Inf]<-0
    annotation_col <- data.frame(
      sample_type = c(rep("group1", each = c_ctrl), rep("group2", each = c_treat))
    )
    rownames(annotation_col) <- colnames(data)
    library("pheatmap")
    pheatmap(
      data,
      #scale = "row",
      color = colorRampPalette(colors = c("#58ACFA", "white", "red"))(10),
      annotation_col = annotation_col,
      cluster_col = F,
      cluster_row = F,
      border = T,
      border_color = "black",
      cellwidth=25,
      cellheight=20,
      annotation_colors = list(sample_type= c(group1 = "red3", group2 = "blue3")),
      filename = fn_stem
    )"""
end


###### Iteratively update housekepping gene
function reoa_update_housekeeping_gene(df_ctrl,
        df_treat,
        names::Vector{String},
        ref_gene::Vector{Bool},
         c_ctrl::Int64,
        c_treat::Int64,
         t_ctrl::Int64,
        t_treat::Int64,
        pval_reo::AbstractFloat,
        padj_deg::AbstractFloat,
         fn_stem::AbstractString,
    ref_gene_new::BitVector, #ref_gene
   iteration_num::Int = 0
	     )
    @time  pvals = compare_reos(
				  Matrix( df_ctrl[!,1:c_ctrl ]),
				  Matrix(df_treat[!,1:c_treat]),
				  Int32(t_ctrl),
				  Int32(t_treat),
				  ref_gene,
                  ref_gene_new)
    delta = pvals[1]
    sd_delta = pvals[2]
    padj = pvals[4]
    pvals_nosingular_matrix_sit =  pvals[5]
    pvals = pvals[3]
    names1 = names[.!ref_gene]
    names1 = names1[pvals_nosingular_matrix_sit]
    pvals_num = (pvals .<= pval_reo)
    padj_num = (padj .<= padj_deg)
    pvals_padj_sit = (pvals_num .& padj_num)
    iteration_1_result = DataFrame()
    insertcols!(iteration_1_result, 1, :names => names1)
    insertcols!(iteration_1_result, 2, :delta => delta)
    insertcols!(iteration_1_result, 3, :sd_delta => sd_delta)
    insertcols!(iteration_1_result, 4, :Pval => pvals)
    insertcols!(iteration_1_result, 5, :Padj => padj)
    iteration_1_fn_out = string(fn_stem, "_iteration_$(iteration_num )_result.tsv")
    CSV.write(iteration_1_fn_out, iteration_1_result, delim = '\t')
    isfile(iteration_1_fn_out) && println("WARN: file $iteration_1_fn_out exists and will be overwritten.")
    print("INFO: N of $(sum(.!pvals_nosingular_matrix_sit)) genes in the $(iteration_num)th calculation is a singular matrix, and there are $(sum(pvals_padj_sit)) differential genes.")
    println("Save the results (gene name, delta, sd_delta, p, FDR) as $(iteration_1_fn_out) file except for the genes of singular matrix.")
    fn_stem_it = string(fn_stem, "_$(iteration_num)")
    delta_p_padj_graph(iteration_1_result,fn_stem_it)
    println("INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the $(iteration_num)th iteration were drawn. Save as ",string(fn_stem_it, "_sd_delta_graph.pdf、"),string(fn_stem_it, "_Pval_graph.pdf、"),string("和",fn_stem_it, "_Padj_graph.pdf"))
    names1 = names1[pvals_padj_sit]
    ref_gene_new_now = (map(x -> ∈(x, names1), names))
    if sum(ref_gene_new_now) == sum(.!ref_gene_new)
        result = DataFrame()
        insertcols!(result, 1, :names => names1)
        insertcols!(result, 2, :delta => delta[pvals_padj_sit])
        insertcols!(result, 3, :sd_delta => sd_delta[pvals_padj_sit])
        insertcols!(result, 4, :Pval => pvals[pvals_padj_sit])
        insertcols!(result, 5, :Padj => padj[pvals_padj_sit])
        fn_out = string(fn_stem, "_result.tsv")
        isfile(fn_out) && println("WARN: file $fn_out exists and will be overwritten.")
        CSV.write(fn_out, result, delim = '\t')
        print("INFO: After $(iteration_num) iterations, the obtained list of differential genes (gene name, delta, sd_delta, p, FDR) is saved as $(fn_out) file, ")
        r_ctrl,  c_ctrl  = size(df_ctrl)
        r_treat, c_treat = size(df_treat)
        deg_exp_ctrl = df_ctrl[((1:r_ctrl)[ref_gene_new_now]),:]
        deg_exp_treat = df_treat[((1:r_treat)[ref_gene_new_now]),:]
        fn_deg_exp = string(fn_stem, "_deg_exp.tsv")
        deg_exp = DataFrame()
        deg_exp[:,:names] = names1
        deg_exp = [deg_exp deg_exp_ctrl deg_exp_treat]
        CSV.write(fn_deg_exp, deg_exp, delim = "\t")
        println("Save differential gene expression profile as $(fn_deg_exp).")
        println("INFO: The final calculation results (gene name, delta, sd_delta, p, FDR) without threshold screening are $(iteration_1_fn_out) files.")
        fn_stem_pheatmap = string(fn_stem, "_deg_exp_graph.pdf")
        deg_exp_graph(deg_exp,fn_stem_pheatmap,c_ctrl,c_treat)
        println("INFO: The expression heat map of DEGs was drawn and saved as ", fn_stem_pheatmap)
        return result
    else
        ref_gene_new = (.!ref_gene_new_now)
        hk_nonhk = copy(names)
        hk_nonhk[ref_gene_new] .= "hk_gene"
        hk_nonhk[.!ref_gene_new] .= "non_hk_gene"
        exp_sit = DataFrame()
        exp_sit[:,:gene_names] = names
        exp_sit[:, :gene_types] = hk_nonhk
        exp_sit = [exp_sit df_ctrl df_treat]
        ref_sum_non = sum((.!ref_gene_new))
        gene_fn_out = string(fn_stem,"_hk_nonhk_gene_",iteration_num,".tsv")
        CSV.write(gene_fn_out, exp_sit, delim = "\t")
        println("INFO: The $(ref_sum_non) non-differential genes obtained in the $(iteration_num)th iteration were added to the list of non-housekeeping genes, and the expression profiles (gene name, housekeeping gene, sample) were saved into $(gene_fn_out).")
        iteration_num = iteration_num + 1
        return reoa_update_housekeeping_gene(df_ctrl,
        df_treat,
           names,
        ref_gene,
          c_ctrl,
         c_treat,
          t_ctrl,
         t_treat,
        pval_reo,
        padj_deg,
         fn_stem,
    ref_gene_new,  #ref_gene
    iteration_num
	     )
    end
end


function reoa(fn_expr::AbstractString = "fn_expr.txt",
        fn_metadata::AbstractString = "fn_metadata.txt";
    expr_threshold::Number = 0,
          pval_reo::AbstractFloat = 0.01,
     pval_sign_reo::AbstractFloat = 1.00,
     padj_sign_reo::AbstractFloat = 0.05,
           hk_file::AbstractString = "$(joinpath(@__DIR__, "..", "hk_gene_file", "HK_genes_info.tsv"))",
           hk_name::AbstractString = "ENSEMBL",
      ref_gene_num::Int = 3000,
    use_housekeeping::AbstractString = "yes",
           species::AbstractString = "human",
    cell_drop_rate::Int = 0,
    gene_drop_rate::Int = 0,
    work_dir::AbstractString = "./",
    use_testdata::AbstractString = "no"
    )
    """
    New implementation.
    """
    check_R_packages()
    cd(work_dir)
    if use_testdata == "yes"
        fn_expr = "$(joinpath(@__DIR__, "..", "test", "fn_expr.txt"))"
        fn_metadata="$(joinpath(@__DIR__, "..", "test", "fn_metadata.txt"))"
    end
    # Import datasets
    isfile(fn_expr) || throw(ArgumentError(fn_expr, " does not exist or is not a regular file."))
    filesize(fn_expr) > 0 || throw(ArgumentError("file for expression matrix has size 0."))
    if species != "human"
        println("INFO: The species information you selected for input data is: ",species)
        if species == "mouse"
            hk_file=replace(hk_file, "HK_genes_info.tsv" => "HK_genes_info_mouse.tsv")
        else
            println("ERROR: Currently, only House Keeping gene information for Human or mouse is provided. Please re-enter Human or mouse.")
        end
    else
        println("INFO: The species information you selected for input data is: ",species)
    end
    fn_stem, = splitext(basename(fn_expr))
    expr = CSV.read(fn_expr, DataFrame)
    a,b=size(expr)
    println("INFO: Successfully read in the expression matrix, including gene names (1st column),  with a size of: ", (a,b))
    if hk_name == "ENTREZID"
        gene_names = [string(expr[i,1]) for i in 1:a]        
    else
        gene_names = convert(Vector{String},expr[!,1])
    end
    expr = expr[!,2:end]
    metadata = CSV.read(fn_metadata, DataFrame)
    md_r,md_c=size(metadata)
    if md_c != 2
        println("ERROR: input ",basename(fn_metadata)," less than 2 columns, please modify and input again.")
    end
    uniq_metadata=unique(metadata[:,2])
    if size(uniq_metadata)[1] != 2
        return(println("ERROR:There are more than two sample groups in the entered group file. Change the sample groups to two groups and enter the sample groups again."))
    end
    group1_samplename=metadata[(1:md_r)[metadata[:,2] .== uniq_metadata[1]],:][:,1]
    group2_samplename=metadata[(1:md_r)[metadata[:,2] .== uniq_metadata[2]],:][:,1]
    df_ctrl  = expr[!, (1:(b-1))[map(x -> ∈(x, group1_samplename), names(expr))] ]
    df_treat = expr[!, (1:(b-1))[map(x -> ∈(x, group2_samplename), names(expr))]]
    no_up_drop_df_ctrl_col = (sum.(eachcol(df_ctrl.>0)) .> cell_drop_rate)
    no_up_drop_df_treat_col = (sum.(eachcol(df_treat.>0)) .> cell_drop_rate)
    println("INFO: There were ",md_r-sum(no_up_drop_df_ctrl_col)-sum(no_up_drop_df_treat_col)," samples with the number of detected genes (non-0 value) less than ",cell_drop_rate,", and the samples were removed.")
    df_ctrl = df_ctrl[:,(1:size(df_ctrl)[2])[no_up_drop_df_ctrl_col]]
    df_treat = df_treat[:,(1:size(df_treat)[2])[no_up_drop_df_treat_col]]
    no_up_drop_df_ctrl_row = (sum.(eachrow(df_ctrl.>0)) .> gene_drop_rate)
    no_up_drop_df_treat_row =(sum.(eachrow(df_treat.>0)) .> gene_drop_rate)
    no_up_drop_row = (no_up_drop_df_ctrl_row .& no_up_drop_df_treat_row)
    df_ctrl = df_ctrl[no_up_drop_row,:]
    df_treat = df_treat[no_up_drop_row,:]
    gene_names = gene_names[no_up_drop_row]
    println("INFO: The number of cells detected in the presence of ",a-sum(no_up_drop_row)," gene (value non-0) was less than ",gene_drop_rate," in the ctrl and/or treat groups, and gene was removed.")
    inds_ctrl = maximum.(eachrow(df_ctrl )) .>= expr_threshold
    inds_treat= maximum.(eachrow(df_treat)) .>= expr_threshold
    inds = inds_ctrl .| inds_treat
    df_ctrl = df_ctrl[inds,:]
    df_treat= df_treat[inds,:]
    gene_names = gene_names[inds]
    println("INFO: Size of the matrices after filtering non-expressed genes")
    r_ctrl,  c_ctrl  = size(df_ctrl)
    r_treat, c_treat = size(df_treat)
    println("INFO: Size of the control group: ",   size(df_ctrl ))
    println("INFO: Size of the treatment group: ", size(df_treat))
    if use_housekeeping == "yes"
        # Read in the list of house-keeping genes
        isfile(hk_file) || throw(ArgumentError(hk_file, " does not exist or is not a regular file."))
        filesize(hk_file) > 0 || throw(ArgumentError("file for house keeping genes has size 0."))
        ref_gene = CSV.read(hk_file, DataFrame)
        ref_gene = ref_gene[!, Symbol(hk_name)]
        c=size(ref_gene)[1]
        if hk_name == "ENTREZID"
            ref_gene = [string(ref_gene[i]) for i in 1:c]
        end
        println("INFO: Successfully read in the house-keeping genes, with a size of: ", (c,1))
        ref_gene = map(x -> ∈(x, ref_gene), gene_names)
        ref_sum_yuan = sum(ref_gene)
        if ref_sum_yuan > ref_gene_num & ref_gene_num < length(ref_gene)
            names_sit = gene_names[ref_gene]
            ref_gene_sit = rand_sit(ref_gene_num,ref_sum_yuan)
            names_sit = names_sit[ref_gene_sit]
            ref_gene = map(x -> ∈(x, names_sit), gene_names)
            ref_sum = length(names_sit)
            print("INFO: The number of housekeeping genes was $(ref_sum_yuan), which was randomly selected and adjusted to $(ref_sum)")
        else
            ref_sum = sum(ref_gene)
            print("The number of housekeeping genes was $(ref_sum_yuan)")
        end
        if ref_sum_yuan == 0
            println(". \nWARN: The expression profile you entered does not contain the housekeeping gene provided by us. Therefore, we do not use housekeeping gene for subsequent analysis.")
            ref_gene = convert(Vector{Bool}, .!(map(x -> ∈(x, gene_names), gene_names)))
            ref_gene_new = (.!ref_gene .| ref_gene)
        else
            hk_nonhk = copy(gene_names)
            hk_nonhk[ref_gene] .= "hk_gene"
            hk_nonhk[.!ref_gene] .= "non_hk_gene"
            exp_sit = DataFrame()
            exp_sit[:,:gene_names] = gene_names
            exp_sit[:, :gene_types] = hk_nonhk
            exp_sit = [exp_sit df_ctrl df_treat]
            ref_sum_non = sum((.!ref_gene))
            gene_fn_out = string(fn_stem,"_hk_nonhk_gene.tsv")
            CSV.write(gene_fn_out, exp_sit, delim = "\t")
            println(", and the number of non-housekeeping genes was $(ref_sum_non). The expression profile (gene name, housekeeping gene or not, and sample) was saved into $(gene_fn_out).")
            c_t_min = min(c_ctrl,c_treat)
            if c_ctrl > 3 .& c_treat > 3
                graph_num = rand_sit(3,c_t_min)
            else
                graph_num = (1:c_t_min)
            end
            for i in graph_num
                hk_non_hk_graph(exp_sit,names(df_ctrl)[i],fn_stem)
                hk_non_hk_graph(exp_sit,names(df_treat)[i],fn_stem)
            end
            ref_gene_new = (ref_gene .& ref_gene)
        end
    else
        ref_gene = convert(Vector{Bool}, .!(map(x -> ∈(x, gene_names), gene_names)))
        ref_gene_new = (.!ref_gene .| ref_gene)
    end
    t_ctrl = get_major_reo_lower_count(c_ctrl , pval_reo)
    t_treat  = get_major_reo_lower_count(c_treat, pval_reo)
    println("INFO: Minimum size for significantly stable REO in the control group: ",   t_ctrl )
    println("INFO: Minimum size for significantly stable REO in the treatment group: ", t_treat)
    result = reoa_update_housekeeping_gene(df_ctrl,
        df_treat,
        gene_names,
        ref_gene,
        c_ctrl,
        c_treat,
        t_ctrl,
        t_treat,
        pval_sign_reo,
        padj_sign_reo,
        fn_stem,
        ref_gene_new
        )
    return result
end

end # module
