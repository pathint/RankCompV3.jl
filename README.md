## RankCompV3 

RankCompV3 is a differential expression analysis algorithm based on relative expression ordering (REO) of gene pairs. It can be applied to bulk or single-cell RNA-sequencing (scRNA-seq) data, microarray gene expression profiles and proteomics profiles, etc. When applied in scRNA-seq data, it can run in single-cell mode or pseudo-bulk mode. The pseudo-bulk mode is expected to improve the accuracy while decreasing runntime and memory cost. 

### Installation

The algorithm is implemented in Julia. Version 1.7 or later is recommended. The simpliest way to install is using the `Pkg` facility in Julia. 

```julia
using Pkg
Pkg.add("RankCompV3")
```

### Examples

#### Quick Start

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
(19999×16 DataFrame
   Row │ Name     pval         padj        n11      n12      n13      n21      n22      n23      n31      n32      n33      Δ1           Δ2          se         z1
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

#### Run your own DEG analysis

You need to prepare two input files before the analysis: metadata file and expression matrix. Both of them should be saved in the `TSV` or `CSV` format and they should be compatible with each other.   

- **metadata file (required).**

 The metadata file contains at least two columns. The first column is the sample names, and the second column is the grouping information. Only two groups are supported at present, therefore, do not include more than two groups. 

 Column names for a metadata should be `Name` and `Group`. 

 See an example metadata file, [fn_metadata.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_metadata.txt).


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

#### Pseudobulk method

For scRNA-seq data, one can carry out a pseudobulk analysis. Rather than using the original single-cell profiles, pseudobulk profiles can be generated and used for DEG analysis. In this method, a random subset of cells from a group is aggregated into a pseudo-bulk profile. 

The pseudobulk method can be turned on by setting `use_pseudobulk = 1`. 

```julia
julia> reoa("scRNA_expr.txt",
    	"scRNA_metadata.txt",
    	use_pseudobulk = 1)
```

By default, the analysis does not use the pseudobulk method (`use_pseudobulk = 0`).  If a value between `6` and `100` is passed to `use_pseudobulk`, that number of pseudobulk profiles will be generated and used for analysis. Other values (except `0` which turns off the pseudobulk method) will generate the default value (`10`) profiles. 


### Optional Parameters

Below lists the optional keyword parameters and their default values.

| Parameter      | Parameter types | Default value       | Parameters to describe                                       |
| -------------- | --------------- | ------------------- | ------------------------------------------------------------ |
| fn_expr        | AbstractString  | "fn_expr.txt"       | Gene expression profile file path. (required)                |
| fn_metadata    | AbstractString  | "fn_metadata.txt"   | Grouping information file path. (required)                   |
| expr_threshold | Number          | 0                   | Gene expression threshold.                                   |
| min_profiles   | Int             | 0                   | Include features (genes) detected in at least this many cells. |
| min_features   | Int             | 0                   | Include profiles (cells) where at least this many features are detected. |
| pval_reo       | AbstractFloat   | 0.01                | Stable threshold for p-value.                                |
| pval_deg       | AbstractFloat   | 0.05                | Significant  reversal threshold for p-value.                 |
| padj_deg       | AbstractFloat   | 0.05                | Significant reversal threshold for FDR  value.               |
| use_pseudobulk | Int             | 0                   | 0 for not using pseudobulk mode, 1 for automatic, 2~5 not used, 6~100 for number of pseudobulk profiles in each group. |
| use_hk_genes   | AbstractString  | "yes"               | Whether to use the housekeeping gene, yes or no.             |
| hk_file        | AbstractString  | "HK_genes_info.tsv" | House-keeping genes  file path.                              |
| gene_name_type | AbstractString  | "ENSEMBL"           | Available choices: Name, REFSEQ, SYMBOL, ENTREZID, ENSEMBL, UNIGENE and GENENAME. |
| ref_gene_max   | Int             | 3000                | If the number of available features is higher than this, take a random sample of this size. |
| ref_gene_min   | Int             | 100                 | If the number is lower than this, ignore the house-keeping genes. |
| n_iter         | Int             | 128                 | Max iterations.                                              |
| n_conv         | Int             | 5                   | Convergence condition: max. difference in the number of DEGs between two consective iterations. |
| work_dir       | AbstractString  | "./"                | Working Directory.                                           |
| use_testdata   | AbstractString  | "no"                | Whether to use the default provided test data for analysis, yes or no. |

#### Example output file: 

- The expression profile of the ctrl group (after preprocessing).  (See [fn_metadata.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_metadata.txt))

- The expression profile of the treat group (after preprocessing).  (See [fn_expr_treat.tsv](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_treat.tsv))

- This file contains Name, pval, padj, n11, n12, n13, n21, n22, n23, n31, n32, n33, Δ1, Δ2, se, z1.  (See [fn_expr_result.tsv](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_result.tsv))

- Graph of Distribution of Expression Values.  (See [fn_expr_expr_dist.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_expr_dist.pdf))

- Heat maps of expression values for the ctrl and treat groups.  (See [fn_expr_expr_heat.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_expr_heat.pdf))

- Distribution of parameters in 3 x 3 contingency tables.  (See [fn_expr_contigency_table.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_contigency_table.pdf))

- Delta distribution.  (See [fn_expr_delta_value.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_delta_value.pdf))

- Distribution of Standard Error (SE).  (See [fn_expr_se.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_se.pdf))

- Distribution of z1.  (See [fn_expr_z1.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_z1.pdf))

- Distribution of p and FDR values.  (See [fn_expr_p_value.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_p_value.pdf))

- Distribution of expression values for DEGs.  (See [fn_expr_degs_expr_dist.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_degs_expr_dist.pdf))

- Heat map of the expression values of DEGs in the ctrl and treat groups.  (See [fn_expr_degs_expr_heat.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_degs_expr_heat.pdf))


##### log file

- [RankCompV3-test-data-output.log](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/RankCompV3_test_data_output.log)
