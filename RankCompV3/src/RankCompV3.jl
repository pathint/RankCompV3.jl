module RankCompV3

using Distributed
using SharedArrays
using Base.Threads
using MultipleTesting
using DataFrames
using DelimitedFiles
using CSV
using Statistics
using RCall
@everywhere using HypothesisTests
@everywhere using Distributions
@everywhere using Random
@everywhere using LinearAlgebra
@everywhere using ArgParse

@everywhere function parsed_args_main(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--fn_expr","-a";help = "Gene expression profile file path";required = true
        "--fn_metadata", "-b";help = "Grouping information file path";required = true
        "--expr_threshold", "-t";arg_type = Int;help = "Gene expression threshold";default = 3
        "--pval_reo", "-p";arg_type = Float64;help = "Stable threshold for p-value";default = 0.01
        "--pval_sign_reo", "-r";arg_type = Float64;help = "Significant reversal threshold for p-value";default = 1.00
        "--padj_sign_reo", "-f";arg_type = Float64;help = "Significant reversal threshold for FDR value";default = 0.05
        "--no_use_housekeeping", "-z";arg_type = Int;help = "Do not use housekeeping gene to set 1, use housekeeping gene to set 0";default = 0
        "--hk_file", "-k";arg_type = String;help = "Housekeeper gene file path";default = "HK_genes_info.tsv"
        "--hk_name", "-e";arg_type = AbstractString;help = "Column name of the column where the housekeeping gene is located";default = "ENSEMBL"
        "--ref_gene_num", "-l";arg_type = Int;help = "The upper limit of the number of housekeeping genes, if it is greater than this value, ref_gene_num housekeeping genes are randomly selected from it";default = 3000
        "--species", "-s";arg_type = AbstractString;help = "Species information of the data, currently supports hunan and mouse";default = "human"
        "--cell_drop_rate", "-c";arg_type = Int;help = "At least how many genes were detected per sample (value non-zero)";default = 0
        "--gene_drop_rate", "-g";arg_type = Int;help = "At least how many cells each gene was detected in (value non-zero)";default = 0
        "--work_dir","-d";arg_type = AbstractString;help = "Working Directory";default = "./"
    end
    return parse_args(s)
end

#   If two genes have the same expression value,
#   a random order is returned.
@everywhere function is_greater(x::Number, y::Number)
	if abs(x - y) < 0.5
		return rand(Bool)    #返回随机bool值，如false或者true
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

@everywhere function sum_reo(
			      c_size::Int32,
			      t_size::Int32,
			      c_sign::Int32,  #稳定对阈值
			      t_sign::Int32,
			     c_count::Int32,  #稳定对数目
			     t_count::Int32
			 )
	"""
	得到三联表
	Return 
	1 = (0, 0，1), a < b in ctrl and a > b in treat
	2 = (0, 1，0), a ~ b in ctrl and a ~ b in treat
    3 = (1, 0，0), a > b in ctrl and a < b in treat
	以此类推
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

@everywhere function compute_pval( c_ctrl::Int32,   # 控制组表达矩阵的列数
    c_treat::Int32,                   # 处理组表达矩阵的列数
     n_ctrl::Int32,                   # ctrl组的显著稳定基因对的阈值
    n_treat::Int32,                   # treat组的显著稳定基因对的阈值
      reo_t::AbstractMatrix{Int32}    # 基因i和管家基因，在ctrl和treat的稳定对数目的两列的表 gene>house 1
   )
    n_ref, = size(reo_t)
    #每个基因和管家基因的显著稳定对的三联表[a b e, c d f, g h i]
    obs_all = reduce(.+, map(x -> sum_reo(c_ctrl, c_treat, n_ctrl, n_treat, x...), eachrow(reo_t)))
    data_3x3 = [obs_all[1]';obs_all[2]';obs_all[3]']
    delta = McCULLACH_test(data_3x3)
    f = open( "delta_nd_z.tsv", "a" )  # 由于一开始可能并没有dat.txt,所以这里用"a"
    m = vcat([delta[1]],[delta[2]])
    writedlm( f,[m] )
	close(f)
    delta = delta[1]
    return delta
end

function compare_reos(ctrl::AbstractMatrix,   # ctrl组的表达矩阵
		     treat::AbstractMatrix,           # treat组的表达矩阵
		    n_ctrl::Int32,                    # Threshold for significantly stable gene pairs   #ctrl组的显著稳定基因对的阈值
		   n_treat::Int32,                    # treat组的显著稳定基因对的阈值
		  ref_gene::Vector{Bool},            # 1-D Bool Array, 0 (false, other gene), 1 (true, ref gene)   #管家基因的0-1列表，0为该基因不是管家基因；1为该基因是管家基因
      ref_gene_new::BitVector
    )
	"""
	Compute p-value for each gene.
	"""
	 r_ctrl,  c_ctrl = size(ctrl)
	r_treat, c_treat = size(treat)
    #判断ctrl组和treat组的基因数是否相相同
	r_ctrl == r_treat || throw(DimensionMismatch("the number of rows of 'ctrl' not equal to the number of rows of 'treat'."))
	#判断ctrl组和管家基因数是否相同
    r_ctrl == length(ref_gene)|| throw(DimensionMismatch("the number of rows of 'ctrl' not equal to the length of 'ref_gene'."))
	#管家基因数
    n_ref = sum(ref_gene)
	#总的样本数
    n_rep = c_ctrl + c_treat
	println("INFO: Number of threads, $(Threads.nthreads())")
    #得到不是管家基因所在的行的位置
	i_deg = (1:r_ctrl)[.!ref_gene]
    #得到管家基因所在的行的位置
	i_ref = (1:r_ctrl)[  (ref_gene_new .|| ref_gene)]
    data = ["delta" "nz_d"]
    writedlm("delta_nd_z.tsv", data, '\t')
    pval_t1 = pmap(i -> compute_pval(Int32(c_ctrl), 
                       Int32(c_treat), 
                       n_ctrl, n_treat, 
                       reduce(vcat, [[Int32(sum(broadcast(is_greater,  ctrl[i,:],
                                        ctrl[j,:])))	Int32(sum(broadcast(is_greater, treat[i,:],
                                        treat[j,:])))] for j in i_ref[(i_ref .!= i)]])),
            #1)对两个分组分别计算稳定对的数目 2）稳定对是对于基因和管家基因 3）遍历每个基因都做一次compute_pval
        #A基因>house keeping gene 1
                i_deg)
    #对奇异矩阵的基因删除
    data = CSV.read("delta_nd_z.tsv", DataFrame)
    pval_t = data[:,1]
    sd_delta = data[:,2]
    pval_t_retain = (pval_t .!= 111)
    pval_t_del = pval_t[pval_t_retain]
    #计算标化后的delta
    sd_delta = sd_delta[pval_t_retain]
    #sd_delta = pval_t_del ./ sqrt(var(pval_t_del))
    #计算p值
    dist = Normal(0,1)
    pval_t_var = pvalue.(dist, sd_delta, tail = :right)
    #对p值计算fdr
    padj = adjust(pval_t_var, BenjaminiHochberg())
    #对奇异矩阵的基因的p和fdr值进行保存，以便后面基因名使用
    pval_t_del = convert(Vector{Float64}, pval_t_del)
	return pval_t_del,sd_delta,pval_t_var,padj,pval_t_retain  #delta,标准化delta,p,padj,非奇异矩阵的位置
end


#输入要为方阵,如#N = [n11 n12 n13, n21 n22 n23, n31 n32 n33]
@everywhere function McCULLACH_test(data::Matrix)
	l_i = size(data)[1]
    l_j = size(data)[2]
    #N,n,r
    N = zeros(l_i-1,l_j-1)
    #D为对角矩阵，n为N的对角线组成的一列矩阵，r为data中仅计算mαβ(同N公式，但仅计算上三角区域)
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
    #行列式,为奇异矩阵，则标记111.0
    if det(N) == 0
        return 111.0,111.0
    else
        N_inv = inv(N)
        #w1,w2,delta1,delta2
        w1 = (D*N_inv*n)/(n'*N_inv*n)
        w2 = N_inv*n
        delta1 = w1'*log.((r .+ 1/2)./(n - r .+ 1/2))
        delta2 = log((1/2 .+ w2' * r)/(1/2 .+ w2' * (n - r)))
        """var_delta1 = 4 * (1 .+ (1/4) * delta1^2)/(n'*N_inv*n)
        var_delta2 = 4 * (1 .+ (1/4) * delta2^2)/(n'*N_inv*n)"""
        var_delta1 = (4 * (1 .+ (1/4) * delta1^2)/(n'*N_inv*n))
        var_delta2 = (4 * (1 .+ (1/4) * delta2^2)/(n'*N_inv*n))
        sd = sqrt((var_delta1 + var_delta2) /2)
        nd_z1 = delta1/sd
        """nd_z2 = delta2/sd
        #计算正态分布的p值
        dist = Normal(0,1)
        p_nd_z1 = pvalue(dist, nd_z1[1], tail = :right)
        p_nd_z2 = pvalue(dist, nd_z2[1], tail = :right)
        #return p_nd_z1"""
        return delta1[1],nd_z1[1]
    end
end

######随机抽取管家基因
@everywhere function rand_sit(ref_gene_num::Int64,ref_sum_yuan::Int64)
    rand_num = unique([rand(1:ref_sum_yuan) for r in 1:ref_gene_num*100])
    if length(rand_num) >= ref_gene_num
        return rand_num[1:ref_gene_num]
    else
        return rand_sit(ref_gene_num)
    end
end

#绘制hk和非hk gene的表达分布图
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

#绘制delta、p和padj的图
function delta_p_padj_graph(data::DataFrame,fn_stem::String)
    @rput data
    @rput fn_stem
    R"""
    library(ggplot2)
    delta_graph=ggplot(data, aes(x=delta)) + 
      geom_histogram(aes(y=..density..),      # 这一步很重要,使用density代替y轴
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
#绘制DEGs的表达热图
function deg_exp_graph(data::DataFrame,fn_stem::String,c_ctrl::Int64,c_treat::Int64)
    @rput data
    @rput fn_stem
    @rput c_ctrl
    @rput c_treat
    R"""
    rownames(data)=data[,1]
    data=data[,-1]
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


######迭代更新housekepping gene
@everywhere function reoa_update_housekeeping_gene(df_ctrl::DataFrame, #ctrl组的表达值谱
        df_treat::DataFrame,      #treat组的表达值谱
        names::Vector{String},  #gene名
        ref_gene::Vector{Bool},   #管家基因0-1bool
         c_ctrl::Int64,           #ctrl组的列数
        c_treat::Int64,           #treat组的列数
         t_ctrl::Int64,   #ctrl查看显著稳定基因对的最低样本数s
        t_treat::Int64,   #treat查看显著稳定基因对的最低样本数s
        pval_reo::AbstractFloat,   # p的阈值
        padj_deg::AbstractFloat,   #fdr值
         fn_stem::AbstractString,  #保存的文件名
    ref_gene_new::BitVector, #ref_gene
   iteration_num::Int = 0  #计算迭代次数
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
    #卡阈值筛选出来的差异基因为1，其他为0
    #卡p
    pvals_num = (pvals .<= pval_reo)
    padj_num = (padj .<= padj_deg)
    pvals_padj_sit = (pvals_num .&& padj_num)
    ###############迭代中的delta############################################
    #保存每次计算的结果，绘制delta分布图
    iteration_1_result = DataFrame()
    insertcols!(iteration_1_result, 1, :names => names1)
    insertcols!(iteration_1_result, 2, :delta => delta)
    insertcols!(iteration_1_result, 3, :sd_delta => sd_delta)
    insertcols!(iteration_1_result, 4, :Pval => pvals)
    # Adjust pvalues for multiple testing
    insertcols!(iteration_1_result, 5, :Padj => padj)
    iteration_1_fn_out = string(fn_stem, "_iteration_$(iteration_num )_result.tsv")
    CSV.write(iteration_1_fn_out, iteration_1_result, delim = '\t')
    isfile(iteration_1_fn_out) && println("WARN: file $iteration_1_fn_out exists and will be overwritten.")
    print("INFO: 第$(iteration_num)次计算中存在$(sum(.!pvals_nosingular_matrix_sit))个基因的N为奇异矩阵,存在$(sum(pvals_padj_sit))个差异基因。")
    println("保留除奇异矩阵的基因外，计算结果（基因名、delta、标化delta、p、fdr）保存为 $(iteration_1_fn_out)文件")
    #绘制delta、p和padj的图
    fn_stem_it = string(fn_stem, "_$(iteration_num)")
    delta_p_padj_graph(iteration_1_result,fn_stem_it)
    println("INFO: 绘制第$(iteration_num)次迭代结果中全基因的sd_delta、Pval和Padj分布图，分别保存为",string(fn_stem_it, "_sd_delta_graph.pdf、"),string(fn_stem_it, "_Pval_graph.pdf、"),string("和",fn_stem_it, "_Padj_graph.pdf"))
    ##################################################################
    #得到差异基因名
    names1 = names1[pvals_padj_sit]
    #差异基因的0，1矩阵，1为是差异基因，0为非差异基因
    ref_gene_new_now = (map(x -> ∈(x, names1), names))
    #if abs(sum(((.!ref_gene_new) .&& ref_gene_new_now)) - sum(.!ref_gene_new)) < 3
    if sum(ref_gene_new_now) == sum(.!ref_gene_new)
        result = DataFrame()
        insertcols!(result, 1, :names => names1)
        insertcols!(result, 2, :delta => delta[pvals_padj_sit])
        insertcols!(result, 3, :sd_delta => sd_delta[pvals_padj_sit])
        insertcols!(result, 4, :Pval => pvals[pvals_padj_sit])
        # Adjust pvalues for multiple testing
        insertcols!(result, 5, :Padj => padj[pvals_padj_sit])
        fn_out = string(fn_stem, "_result.tsv")
        isfile(fn_out) && println("WARN: file $fn_out exists and will be overwritten.")
        CSV.write(fn_out, result, delim = '\t')
        print("INFO: 经过 $(iteration_num) 次迭代，得到的差异基因列表（基因名、delta、标化delta、p、fdr）保存为 $(fn_out) 文件，")
        #输出差异基因表达，exp_deg
        r_ctrl,  c_ctrl  = size(df_ctrl)
        r_treat, c_treat = size(df_treat)
        deg_exp_ctrl = df_ctrl[((1:r_ctrl)[ref_gene_new_now]),:]
        deg_exp_treat = df_treat[((1:r_treat)[ref_gene_new_now]),:]
        fn_deg_exp = string(fn_stem, "_deg_exp.tsv")
        deg_exp = DataFrame()
        deg_exp[:,:names] = names1
        deg_exp = [deg_exp deg_exp_ctrl deg_exp_treat]
        CSV.write(fn_deg_exp, deg_exp, delim = "\t")
        println("差异基因表达谱保存为：$(fn_deg_exp)")
        println("INFO: 不做阈值筛选的最终计算结果（基因名、delta、标化delta、p、fdr）为 $(iteration_1_fn_out)文件")
        fn_stem_pheatmap = string(fn_stem, "_deg_exp_graph.pdf")
        deg_exp_graph(deg_exp,fn_stem_pheatmap,c_ctrl,c_treat)
        println("INFO: 绘制DEGs的表达热图，保存为", fn_stem_pheatmap)
        return result
    else
        ref_gene_new = (.!ref_gene_new_now)
        #############################################看下管家基因和非管家基因表达谱
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
        println("INFO: 第$(iteration_num)次迭代的非管家基因数为$(ref_sum_non)，表达谱(基因名、是否管家基因、样本）已存入为$(gene_fn_out)")
        #####################################################################################
        iteration_num = iteration_num + 1
        return reoa_update_housekeeping_gene(df_ctrl,
        df_treat,
           names,
        ref_gene,
          c_ctrl,
         c_treat,
          t_ctrl,  #ctrl查看显著稳定基因对的最低样本数s
         t_treat,  #treat查看显著稳定基因对的最低样本数s
        pval_reo,  #基因的p值
        padj_deg,  #fdr值
         fn_stem,  #文件名带有exp的名字开头
    ref_gene_new,  #ref_gene
    iteration_num  #迭代次数
	     )
    end
end
function reoa(fn_expr::AbstractString="expr.txt",        # 基因表达谱文件路径
        fn_metadata::AbstractString="metadata.txt";      # 分组信息文件路径
    expr_threshold::Number = 3,           # 基因表达阈值
       #n_subsample::Int = 100,            # 如果分组中样本数大于该数值，则随机抽取n_subsample个样本
          pval_reo::AbstractFloat = 0.01,   # 稳定对p值的阈值
     pval_sign_reo::AbstractFloat = 1.00,   # 显著逆转对p值的阈值
     padj_sign_reo::AbstractFloat = 0.05,   # 显著逆转对FDR值的阈值
           hk_file::AbstractString = "HK_genes_info.tsv",  # 管家基因文件路径
           hk_name::AbstractString = "ENSEMBL",            # 管家基因所在列的列名
      ref_gene_num::Int = 3000,             # 管家基因数目上限，如果大于该值.则从中随机抽取ref_gene_num个管家基因
    no_use_housekeeping::Int = 0,           # 不使用管家基因设置1，使用管家基因设置0
           species::AbstractString = "human",   # 数据的物种信息，目前支持hunan和mouse
    cell_drop_rate::Int = 0,                # 每个样本至少有多少个基因被检测到（值非零）
    gene_drop_rate::Int = 0,               # 每个基因至少在多少个细胞中被检测到（值非零）
    work_dir::AbstractString = "./"         #设置工作目录
    )
    """
    New implementation.
    """
    cd(work_dir)
    # Import datasets
    isfile(fn_expr) || throw(ArgumentError(fn_expr, " does not exist or is not a regular file."))
    filesize(fn_expr) > 0 || throw(ArgumentError("file for expression matrix has size 0."))
    if species != "human"
        println("INFO: 你选择输入数据的物种信息为：",species)
        if species == "mouse"
            hk_file=replace(hk_file, "HK_genes_info.tsv" => "HK_genes_info_mouse.tsv")
        else
            println("ERROR: 目前仅提供human或mouse的house keeping gene信息，请重新输入human或mouse")
        end
    else
        println("INFO: 你选择输入数据的物种信息为：",species)
    end
    fn_stem, = splitext(basename(fn_expr))
    expr = CSV.read(fn_expr, DataFrame)
    a,b=size(expr)
    println("INFO: Successfully read in the expression matrix, including gene names (1st column),  with a size of: ", (a,b))
    # Drop rows with missing values (NA), dropmissing!(df)
    # or replace zeros. 
    #dropmissing!(expr, disallowmissing=true)
    # assuming the first column are the gene names
    if hk_name == "ENTREZID"
        gene_names = [string(expr[i,1]) for i in 1:a]        
    else
        gene_names = convert(Vector{String},expr[!,1])
    end
    """if length(unique(gene_names)) != length(gene_names)
        return print("你输入的表达谱中存在相同的基因$hk_name")
    end"""
    expr = expr[!,2:end]
    metadata = CSV.read(fn_metadata, DataFrame)
    uniq_metadata=unique(metadata[:,2])
    if size(uniq_metadata)[1] != 2
        return(println("ERROR:输入的样本分组大于两组，请修改为两个分组后重新输入分组文件。"))
    end
    #提取两个分组的样本信息
    group1_samplename=metadata[(1:(b-1))[metadata[:,2] .== uniq_metadata[1]],:][:,1]
    group2_samplename=metadata[(1:(b-1))[metadata[:,2] .== uniq_metadata[2]],:][:,1]
    # Filter, filter(row -> maximum(row)>5, df)
    #在表达谱中对分组样本信息进行分类样本
    df_ctrl  = expr[!, (1:(b-1))[map(x -> ∈(x, group1_samplename), names(expr))] ]
    df_treat = expr[!, (1:(b-1))[map(x -> ∈(x, group2_samplename), names(expr))]]
    #每个样本至少有cell_drop_rate个基因被检测到（值非零）
    no_up_drop_df_ctrl_col = (sum.(eachcol(df_ctrl.>0)) .> cell_drop_rate)
    no_up_drop_df_treat_col = (sum.(eachcol(df_treat.>0)) .> cell_drop_rate)
    println("INFO: 存在",b-sum(no_up_drop_df_ctrl_col)-sum(no_up_drop_df_treat_col),"个样本检测到的基因数（值非0）少于",cell_drop_rate,"，并对样本进行去除。")
    df_ctrl = df_ctrl[:,(1:size(df_ctrl)[2])[no_up_drop_df_ctrl_col]]
    df_treat = df_treat[:,(1:size(df_treat)[2])[no_up_drop_df_treat_col]]
    #每个基因至少在gene_drop_rate个细胞中被检测到（值非零）
    no_up_drop_df_ctrl_row = (sum.(eachrow(df_ctrl.>0)) .> gene_drop_rate)
    no_up_drop_df_treat_row =(sum.(eachrow(df_treat.>0)) .> gene_drop_rate)
    no_up_drop_row = (no_up_drop_df_ctrl_row .&& no_up_drop_df_treat_row)
    df_ctrl = df_ctrl[no_up_drop_row,:]
    df_treat = df_treat[no_up_drop_row,:]
    gene_names = gene_names[no_up_drop_row]
    println("INFO: 存在",a-sum(no_up_drop_row),"个gene检测到的细胞数（值非0）在ctrl和/或treat组中少于",gene_drop_rate,"，并对gene进行去除。")
    #依次判断df_ctrl是否大于expr_threshold。对表达值取阈值
    inds_ctrl = maximum.(eachrow(df_ctrl )) .>= expr_threshold
    inds_treat= maximum.(eachrow(df_treat)) .>= expr_threshold
    inds = inds_ctrl .| inds_treat
    df_ctrl = df_ctrl[inds,:]
    df_treat= df_treat[inds,:]
    gene_names = gene_names[inds]
    #println("Size of the matrix after filtering non-expressed genes: ", size(expr))
    println("INFO: Size of the matrices after filtering non-expressed genes")
    r_ctrl,  c_ctrl  = size(df_ctrl)
    r_treat, c_treat = size(df_treat)
    println("INFO: Size of the control group: ",   size(df_ctrl ))
    println("INFO: Size of the treatment group: ", size(df_treat))
    #判断是否使用hk基因
    if no_use_housekeeping == 0
        # Read in the list of house-keeping genes
        isfile(hk_file) || throw(ArgumentError(hk_file, " does not exist or is not a regular file."))
        filesize(hk_file) > 0 || throw(ArgumentError("file for house keeping genes has size 0."))
        #读入管家基因列表
        ref_gene = CSV.read(hk_file, DataFrame)
        #提取出hk_name列名的列，默认为"ENSEMBL"
        ref_gene = ref_gene[!, Symbol(hk_name)]
        c=size(ref_gene)[1]
        if hk_name == "ENTREZID"
            ref_gene = [string(ref_gene[i]) for i in 1:c]
        end
        println("INFO: Successfully read in the house-keeping genes, with a size of: ", (c,1))
        ref_gene = map(x -> ∈(x, ref_gene), gene_names)
        ################对管家基因数卡阈值
        ref_sum_yuan = sum(ref_gene)
        if ref_sum_yuan > ref_gene_num && ref_gene_num < length(ref_gene)
            #得到表达谱的管家基因
            names_sit = gene_names[ref_gene]
            #随机抽取指定数目的管家基因
            ref_gene_sit = rand_sit(ref_gene_num,ref_sum_yuan)
            names_sit = names_sit[ref_gene_sit]
            ref_gene = map(x -> ∈(x, names_sit), gene_names)
            ref_sum = length(names_sit)
            print("INFO: 总管家基因数目为$(ref_sum_yuan)个，随机抽取后调整为$(ref_sum)个,")
        else
            ref_sum = sum(ref_gene)
            print("总管家基因数目为$(ref_sum_yuan)个,")
        end
        #保存管家基因和非管家基因表达谱
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
        println("非管家基因数为$(ref_sum_non)，表达谱(基因名、是否管家基因、样本）已存入为$(gene_fn_out)")
        #随机抽取样本，绘制管家基因和非管家基因表达图
        c_t_min = min(c_ctrl,c_treat)
        if c_ctrl > 3 .&& c_treat > 3
            graph_num = rand_sit(3,c_t_min)
        else
            graph_num = (1:c_t_min)
        end
        for i in graph_num
            hk_non_hk_graph(exp_sit,names(df_ctrl)[i],fn_stem)
            hk_non_hk_graph(exp_sit,names(df_treat)[i],fn_stem)
        end
        #############################
        ref_gene_new = (ref_gene .&& ref_gene)
    else
        ref_gene = convert(Vector{Bool}, .!(map(x -> ∈(x, gene_names), gene_names)))
        ref_gene_new = (.!ref_gene .|| ref_gene)
    end
    #############################
    """# 100 is arbitrary
    sum(ref_gene) > 100 || throw(ErrorException("Too few house-keeping genes in datasets. It must be > 100."))
    # Using Bionomial distribution to obtain threshold value 
    sub_ctrl = 0
    sub_treat= 0
    #对于大样本简化样本量，简化计算
    if c_ctrl <= n_subsample
        #get_major_reo_lower_count为查看显著稳定基因对的最低样本数s
        t_ctrl = get_major_reo_lower_count(c_ctrl , pval_reo)
    else
        println("INFO: Subsampling is turned on for the control group.")
        t_ctrl = get_major_reo_lower_count(n_subsample, pval_reo)
        sub_ctrl = n_subsample
    end
    if c_treat <= n_subsample
        t_treat  = get_major_reo_lower_count(c_treat, pval_reo)
    else
        println("INFO: Subsampling is turned on for the treatment group.")
        t_treat  = get_major_reo_lower_count(n_subsample, pval_reo)
        sub_treat= n_subsample
    end"""
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
        t_ctrl,   #ctrl查看显著稳定基因对的最低样本数s
        t_treat,  #treat查看显著稳定基因对的最低样本数s
        pval_sign_reo, #基因的p值
        padj_sign_reo, #fdr值
        fn_stem,   #文件名
        ref_gene_new
        )
    return result
end

parsed_args = parsed_args_main(ARGS)
function julia_main()
    reoa(parsed_args["fn_expr"],                            # 基因表达谱文件路径
            parsed_args["fn_metadata"];                 # 分组信息文件路径
        expr_threshold = parsed_args["expr_threshold"],     # 基因表达阈值
            pval_reo = parsed_args["pval_reo"],      # 稳定对p的阈值
        pval_sign_reo = parsed_args["pval_sign_reo"], # 显著逆转对p的阈值
        padj_sign_reo = parsed_args["padj_sign_reo"], # 显著逆转对FDR的阈值
            hk_file = parsed_args["hk_file"],      # 管家基因文件
            hk_name = parsed_args["hk_name"],      # 管家基因所在列的列名
        ref_gene_num = parsed_args["ref_gene_num"],            # 管家基因数目上限，如果大于.则从中随机抽取ref_gene_num个管家基因
        no_use_housekeeping = parsed_args["no_use_housekeeping"],# 不使用管家基因设置1，使用管家基因设置0
            species = parsed_args["species"],      # 数据的物种信息，目前支持hunan和mouse
        cell_drop_rate = parsed_args["cell_drop_rate"],          # 每个样本至少有多少个基因被检测到（值非零）
        gene_drop_rate = parsed_args["gene_drop_rate"],          # 每个基因至少在多少个细胞中被检测到（值非零）
        work_dir=parsed_args["work_dir"]                         #工作目录
        )
    return 0 # if things finished successfully
end


end # module
