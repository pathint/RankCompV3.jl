# RankCompV3 

RankCompV3 is a differential expression analysis algorithm based on relative expression ordering (REO) of gene pairs. It can be applied to bulk or single-cell RNA-sequencing (scRNA-seq) data, microarray gene expression profiles and proteomics profiles, etc. When applied in scRNA-seq data, it can run in single-cell mode or pseudo-bulk mode. The pseudo-bulk mode is expected to improve the accuracy while decreasing runntime and memory cost. 

## 1 Used in the Julia language

### 1.1 Installation

The algorithm is implemented in Julia. Version 1.7 or later is recommended. The simpliest way to install is using the `Pkg` facility in Julia. 

```julia
using Pkg
Pkg.add("RankCompV3")
```

### 1.2 Examples

#### 1.2.1 Quick Start

Run a test job with the input files distributed with the package.

```julia
julia> using RankCompV3
# Use the default values for the following other parameters. If you need to modify the parameters, add them directly.
julia> result = reoa(use_testdata="yes")
```

The analysis results and a few plots will be generated and saved in the current work directory. They are also returned by the `reoa` function and can be captured by assign the returned values to a variable,  e.g., `result` in the above example.  

The first return value is a DataFrame, where rows are genes and columns are statistical values for each gene. All the genes passing the basic preprocessing step are retained. 


```julia
julia> result
19999×2 DataFrame
   Row │ gene_name  group1_vs_group2
       │ Any        Any
───────┼─────────────────────────────
     1 │ DE1        up
     2 │ DE2        up
     3 │ DE3        up
     4 │ DE4        up
     5 │ DE5        up
     6 │ DE6        up
   ⋮   │     ⋮             ⋮
 19995 │ EE19996    up
 19996 │ EE19997    up
 19997 │ EE19998    up
 19998 │ EE19999    up
 19999 │ EE20000    up
                   19988 rows omitted
```

#### 1.2.2 Run your own DEG analysis

You need to prepare two input files before the analysis: metadata file and expression matrix. `.rds`, `.csv`, `.txt`, `.tsv` and `.RData` files are supported. 

- **metadata file (required).**

 The metadata file contains at least two columns. The first column is the sample names, and the second column is the grouping information. Only two groups are supported at present, therefore, do not include more than two groups. 

 Column names for a metadata should be `Name` and `Group`. 

 See an example metadata file, [fn_meta.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_meta.txt).


- **expression matrix file (required).**

 The first column is the gene name and the column header should be `Name` and the rest columns are profiles for each cell or each sample. Each column header should be the sample name which appears in the metadata file.

 See an example expression matrix file, [fn_expr.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_expr.txt).

 Raw counts or expression values are recommended to use. Other values, e.g, FPKM, RPKM, TPM, log(counts) and log(normalized counts), can also be used, though normalization and batch effect removal are neither necessary nor recommended. 

Once the files are ready, you can carry out the DEG analysis with the default settings as follows. 

```julia
julia> using RankCompV3
# Use the default values for the following other parameters. If you want to modify the parameters, add them directly.
julia> reoa("/public/yanj/data/fn_expr.txt",
	"/public/yanj/data/fn_metadata.txt")
```

Other parameters can be set by passing the value to the corresponding keyword. 

```julia
julia> reoa("/public/yanj/data/fn_expr.txt",
    	"/public/yanj/data/fn_meta.txt",
    	expr_threshold = 0,
    	min_profiles = 0,
    	min_features = 0,
    	pval_reo = 0.01,
     	pval_deg = 0.05,
     	padj_deg = 0.05,
    	n_pseudo = 0,
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

#### 1.2.3 Pseudobulk method

For scRNA-seq data, one can carry out a pseudobulk analysis. Rather than using the original single-cell profiles, pseudobulk profiles can be generated and used for DEG analysis. In this method, a random subset of cells from a group is aggregated into a pseudo-bulk profile. 

The pseudobulk method can be turned on by setting `n_pseudo > 0`. 

```julia
julia> reoa("scRNA_expr.txt",
    	"scRNA_metadata.txt",
    	n_pseudo = 1)
```

By default, profiling does not use the pseudo-bulk method (`n_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate the number of samples in each group after pseudo-bulk is combined.


### 1.3 Optional Parameters

Below lists the optional keyword parameters and their default values.

| Parameter      | Parameter types | Default value       | Parameters to describe                                       |
| -------------- | --------------- | ------------------- | ------------------------------------------------------------ |
| fn_expr        | AbstractString  | "fn_expr.txt"       | Gene expression profile file path. (required).               |
| fn_metadata    | AbstractString  | "fn_metadata.txt"   | Grouping information file path. (required)                   |
| min_profiles   | Int             | 0                   | Include features (genes) detected in at least this many cells. |
| min_features   | Int             | 0                   | Include profiles (cells) where at least this many features are detected. |
| pval_reo       | AbstractFloat   | 0.01                | Stable threshold for p-value.                                |
| pval_deg       | AbstractFloat   | 0.05                | Significant  reversal threshold for p-value.                 |
| padj_deg       | AbstractFloat   | 0.05                | Significant reversal threshold for FDR  value.               |
| n_pseudo       | Int             | 0                   | 0 indicates that the pseudo-bulk mode is not used, and other values indicate the number of samples in each group after pseudo-bulk is combined. |
| use_hk_genes   | AbstractString  | "yes"               | Whether to use the housekeeping gene, yes or no.             |
| hk_file        | AbstractString  | "HK_genes_info.tsv" | House-keeping genes  file path.                              |
| gene_name_type | AbstractString  | "ENSEMBL"           | Available choices: Name, REFSEQ, SYMBOL, ENTREZID, ENSEMBL, UNIGENE and GENENAME. |
| ref_gene_max   | Int             | 3000                | If the number of available features is higher than this, take a random sample of this size. |
| ref_gene_min   | Int             | 100                 | If the number is lower than this, ignore the house-keeping genes. |
| n_iter         | Int             | 128                 | Max iterations.                                              |
| n_conv         | Int             | 5                   | Convergence condition: max. difference in the number of DEGs between two consective iterations. |
| work_dir       | AbstractString  | "./"                | Working Directory.                                           |
| use_testdata   | AbstractString  | "no"                | Whether to use the default provided test data for analysis, yes or no. |

### 1.4 Example output file

#### 1.4.1 result

- The expression profile (after preprocessing).  (See [fn_expr_df_expr.tsv](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_df_expr.tsv))
- The meta data (after preprocessing).  (See [fn_expr_df_meta.tsv](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_df_meta.tsv))
- This file contains Name, pval, padj, n11, n12, n13, n21, n22, n23, n31, n32, n33, Δ1, Δ2, se, z1, up_down.  (See [fn_expr_group1_group2_result.tsv](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_result.tsv))
- Up-down-regulation of all genes in all groups (up is up-regulation, down is down, no change means that the gene is not recognized as up-regulation/down-regulation) (See [fn_expr_gene_up_down.tsv](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_gene_up_down.tsv))
- Graph of Distribution of Expression Values.  (See [fn_expr_group1_group2_expr_dist.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_expr_dist.pdf))
- Heat maps of expression values for the ctrl and treat groups.  (See [fn_expr_group1_group2_expr_heat.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_expr_heat.pdf))
- Distribution of parameters in 3 x 3 contingency tables.  (See [fn_expr_group1_group2_contigency_table.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_contigency_table.pdf))
- Delta distribution.  (See [fn_expr_group1_group2_delta_value.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_delta_value.pdf))
- Distribution of Standard Error (SE).  (See [fn_expr_group1_group2_se.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_se.pdf))
- Distribution of z1.  (See [fn_expr_group1_group2_z1.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_z1.pdf))
- Distribution of p and FDR values.  (See [fn_expr_group1_group2_p_value.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_p_value.pdf))
- Distribution of expression values for DEGs.  (See [fn_expr_group1_group2_degs_expr_dist.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_degs_expr_dist.pdf))
- Heat map of the expression values of DEGs in the ctrl and treat groups.  (See [fn_expr_group1_group2_degs_expr_heat.pdf](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/fn_expr_group1_group2_degs_expr_heat.pdf))


#### 1.4.2 log file

- [RankCompV3_test_data_output.log](https://github.com/yanjer/testdata-output/blob/master/RankCompV3-test-data-output/RankCompV3_test_data_output.log)

## 2 Used in the R language

### 2.1 Installation

##### 2.1.1 You can install just like any other R packages by `JuliaCall`

```R
install.packages("JuliaCall")
```

##### 2.1.2 To use you must have a working installation of Julia. This can be easily done via: `JuliaCall`

```R
library(JuliaCall)
install_julia()
```

##### 2.1.3 which will automatically install and setup a version of Julia specifically for use with `JuliaCall`. Or you can do

```R
library(JuliaCall)
julia <-julia_setup()
```

##### 2.1.4 Download RankCompV3

```julia
julia_install_package_if_needed("RankCompV3")
```

### 2.2 Examples

#### 2.2.1 Quick Start

Run a test job with the input files distributed with the package.

```R
julia_library("RankCompV3")
result <- julia_do.call("reoa",list(use_testdata="yes"),need_return="Julia",show_value=FALSE)
```

The analysis results and a few plots will be generated and saved in the current work directory. They are also returned by the `reoa` function and can be captured by assign the returned values to a variable,  e.g., `result` in the above example.  

The first return value is a DataFrame, where rows are genes and columns are statistical values for each gene. All the genes passing the basic preprocessing step are retained. 


```R
> result
Julia Object of type Tuple{DataFrames.DataFrame, DataFrames.DataFrame, DataFrames.DataFrame}.
19999×2 DataFrame
   Row │ gene_name  group1_vs_group2
       │ Any        Any
───────┼─────────────────────────────
     1 │ DE1        up
     2 │ DE2        up
     3 │ DE3        up
     4 │ DE4        up
     5 │ DE5        up
     6 │ DE6        up
   ⋮   │     ⋮             ⋮
 19995 │ EE19996    up
 19996 │ EE19997    up
 19997 │ EE19998    up
 19998 │ EE19999    up
 19999 │ EE20000    up
                   19988 rows omitted
```

#### 2.2.2 Run your own DEG analysis

You need to prepare two input files before the analysis: metadata file and expression matrix. You need to prepare two input files before the analysis: metadata file and expression matrix. `.rds`, `.csv`, `.txt`, `.tsv` and `.RData` files are supported.    

- **metadata file (required).**

 The metadata file contains at least two columns. The first column is the sample names, and the second column is the grouping information. Only two groups are supported at present, therefore, do not include more than two groups. 

 Column names for a metadata should be `Name` and `Group`. 

 See an example metadata file, [fn_metadata.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_metadata.txt).


- **expression matrix file (required).**

 The first column is the gene name and the column header should be `Name` and the rest columns are profiles for each cell or each sample. Each column header should be the sample name which appears in the metadata file.

 See an example expression matrix file, [fn_expr.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_expr.txt).

 Raw counts or expression values are recommended to use. Other values, e.g, FPKM, RPKM, TPM, log(counts) and log(normalized counts), can also be used, though normalization and batch effect removal are neither necessary nor recommended. 

Once the files are ready, you can carry out the DEG analysis with the default settings as follows. 

```R
julia_library("RankCompV3")
# Use the default values for the following other parameters. If you want to modify the parameters, add them directly.
julia_do.call("reoa",list("/public/yanj/data/fn_expr.txt",
						"/public/yanj/data/fn_meta.txt"),need_return="Julia",show_value=FALSE)
```

Other parameters can be set by passing the value to the corresponding keyword. 

```R
julia_do.call("reoa",list("/public/yanj/data/fn_expr.txt",
    	"/public/yanj/data/fn_meta.txt",
    	expr_threshold = 0,
    	min_profiles = 0,
    	min_features = 0,
    	pval_reo = 0.01,
     	pval_deg = 0.05,
     	padj_deg = 0.05,
    	n_pseudo = 0,
    	use_hk_genes = "yes"
    	hk_file = "HK_genes_info.tsv",
    	gene_name_type = "ENSEMBL",
    	ref_gene_max = 3000,
    	ref_gene_min = 100
    	n_iter = 128,
    	n_conv = 5,
    	work_dir = "./",
    	use_testdata = "no"),need_return="Julia",show_value=FALSE)
```

#### 2.2.3 Pseudobulk method

For scRNA-seq data, one can carry out a pseudobulk analysis. Rather than using the original single-cell profiles, pseudobulk profiles can be generated and used for DEG analysis. In this method, a random subset of cells from a group is aggregated into a pseudo-bulk profile. 

The pseudobulk method can be turned on by setting `n_pseudo = 1`. 

```R
julia_do.call("reoa",list("scRNA_expr.txt",
						"scRNA_metadata.txt",
        				n_pseudo = 1),need_return="Julia",show_value=FALSE)
```

By default, profiling does not use the pseudo-bulk method (`n_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate the number of samples in each group after pseudo-bulk is combined.


### 2.3 Optional Parameters

See 1.3 Optional Parameters.

### 2.4 Example output file

See 1.4 Example output file.

## 3 Specific experimental designs

In this chapter, we will present some examples of the use of experimental designs.

### 3.1 Datasets from different sources

Differential expression analysis often involves datasets from different sources, such as different sequencing platforms, different sample sources, or different sample processing methods. RankCompV3 supports direct column concatenation of expression profiles for integrated analysis of different datasets.

```julia
julia> expr
5×11 DataFrame
 Row │ gene_name  s1     s2     s3     c1     c8     c6     a1     a2     b1     b6
     │ String7    Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64
─────┼─────────────────────────────────────────────────────────────────────────────────
   1 │ DE1           10     53     18     84     18    200    180    176     26     24
   2 │ DE2            5     22     59     11     33     26     37     36     40     79
   3 │ DE3           62     39     18     19      8    157     81     89    220     97
   4 │ DE4            6    116    131      3     49     15    205     63     75    228
   5 │ DE5           28     31     27     58     16    505     92    119    157    261
```

Among them, s, c, a and b samples are data from different sources. Samples s and c are in the same group, and samples a and b are in the same group.

```julia
julia> meta
10×2 DataFrame
 Row │ sample_name  group
     │ String15     String7
─────┼──────────────────────
   1 │ s1           group1
   2 │ s2           group1
   3 │ s3           group1
   4 │ c1           group1
   5 │ c8           group1
   6 │ c6           group1
   7 │ a1           group2
   8 │ a2           group2
   9 │ b1           group2
  10 │ b6           group2
```

### 3.2 Two or more groups

When multiple cell types are included, the characteristics of each cell class need to be analyzed and each cell class compared to the rest. Suppose there are three types of cells, such as cell types s, c, a, and b.

```julia
julia> meta
10×2 DataFrame
 Row │ sample_name  group
     │ String       String
─────┼────────────────────────
   1 │ s1           celltype1
   2 │ s2           celltype1
   3 │ s3           celltype1
   4 │ c1           celltype2
   5 │ c8           celltype2
   6 │ c6           celltype2
   7 │ a1           celltype3
   8 │ a2           celltype3
   9 │ b1           celltype4
  10 │ b6           celltype4
```

By default, conditions are listed in order of precedence.

```julia
julia> unique(meta.group)
4-element Vector{String}:
 "celltype1"
 "celltype2"
 "celltype3"
 "celltype4"
```

### 3.3 Combined sample

For analyzing single-cell data, drop-out phenomenon often exists. The `n_pseudo` parameter of pseudobulk mode of RankCompV3 can be used to combine different cells to reduce the effect of drop-out.

```julia
julia> reoa("/public/yanj/data/fn_expr.txt",
    	"/public/yanj/data/fn_meta.txt",
    	n_pseudo = 10)
```

By default, profiling does not use the pseudo-bulk method (`n_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate the number of samples in each group after pseudo-bulk is combined.























