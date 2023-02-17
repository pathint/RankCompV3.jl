## RankCompV3

​	RankCompV3 is a differentially expressed gene recognition algorithm based on relative expression order relation REO. The tool is developed based on the julia language. The details are described below. julia recommends using version 1.7 or later.

**RankCompV3 package in julia：https://github.com/pathint/RankCompV3.jl**

### Scope of application of RankCompV3

| Classification                              | Applicable scope                                             |
| ------------------------------------------- | ------------------------------------------------------------ |
| **Species**                                 | Unrestricted species                                         |
| **Species supported by housekeeping genes** | human                                                        |
| **Data type**                               | Microarray data, high-throughput sequencing data (RNA-seq), single-cell sequencing data (scRNA-seq), Cell expression by linear amplification and sequencing (CEL-seq) , methylation data and proteome data, etc |
| **Data format**                             | Count/log2Count/FPKM/RPKM/TPM/Expression values, etc         |
| **Genotype**                                | REFSEQ/SYMBOL/ENTREZID/ENSEMBL/UNIGENE/GENENAME              |

### Tool Dependency package

​	**Where julia related software packages have been added to the environment. **

```julia
CSV v0.10.9
CategoricalArrays v0.10.7
DataFrames v1.4.4
Distributions v0.25.80
HypothesisTests v0.10.11
MultipleTesting v0.5.1
Plots v1.38.4
StatsBase v0.33.21
StatsPlots v0.15.4
DelimitedFiles
Distributed
LinearAlgebra
Random
SharedArrays
Statistics
```

### Configure

#### RankCompV3 package in julia

##### Use directly

```shell
#Load Pkg in julia
using Pkg
#The RankCompV3 package is required for the first use
Pkg.add("RankCompV3")
#or
Pkg.add(url="https://github.com/pathint/RankCompV3.jl.git")
```

##### Local run usage

```shell
#configured in linux
#clone the RankCompV3 package from github to local
git clone https://github.com/pathint/RankCompV3.jl.git
#Load the project dependency package
#[path] is the path from git clone https://github.com/yanjer/RankCompV3.git
julia --project=RankCompV3
#Download dependency packages
##Enter Pkg mode. ]
instantiate
```

### Parameters for details

| Parameter      | Parameter types              | Default value     | Parameters to describe                                       |
| -------------- | ---------------------------- | ----------------- | ------------------------------------------------------------ |
| fn_expr        | AbstractString               | fn_expr.txt       | Gene expression profile file path. (required)                |
| fn_metadata    | AbstractStringAbstractString | fn_metadata.txt   | Grouping information file path. (required)                   |
| expr_threshold | NumberNumber                 | 0                 | Gene expression threshold.                                   |
| min_profiles   | Int                          | 0                 | Include features (genes) detected in at least this many cells. |
| min_features   | Int                          | 0                 | Include profiles (cells) where at least this many features are detected. |
| pval_reo       | AbstractFloatAbstractFloat   | 0.01              | Stable threshold for p-value.                                |
| pval_deg       | AbstractFloatAbstractFloat   | 0.05              | Significant  reversal threshold for p-value.                 |
| padj_deg       | AbstractFloatAbstractFloat   | 0.05              | Significant reversal threshold for FDR  value.               |
| use_pseudobulk | Int                          | 0                 | 0 for not using pseudobulk mode, 1 for automatic, 2~5 not used, 6~100 for number of pseudobulk profiles in each group. |
| use_hk_genes   | AbstractString               | yes               | Whether to use the housekeeping gene, yes or no.             |
| hk_file        | AbstractString               | HK_genes_info.tsv | Housekeeper gene  file path.                                 |
| gene_name_type | AbstractString               | ENSEMBL           | Available choices: ENSEMBL, Symbol, ENTREZID ...             |
| ref_gene_max   | Int                          | 3000              | If the number of available features is higher than this, take a random sample of this size. |
| refinement     | Int                          | 100               | If the number is lower than this, ignore the house-keeping genes. |
| n_iter         | Int                          | 128               | Max iterations.                                              |
| n_conv         | Int                          | 5                 | Convergence condition: max. difference in the number of DEGs between two consective iterations |
| work_dir       | AbstractString               | ./                | Working Directory.                                           |
| use_testdata   | AbstractString               | no                | Whether to use the default provided test data for analysis, yes or no. |

### Input File Description

| File                        | Description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| **Expression profile file** | Each row represents the gene, each column represents the sample and the expression matrix of the gene in the first column. |
| **Metadata file**           | First column sample name, second column group information.   |
| **housekeeping gene file**  | Housekeeping gene file with built-in support for several gene forms: Name, REFSEQ, SYMBOL, ENTREZID, ENSEMBL, UNIGENE and GENENAME. Also supports external incoming or not using housekeeping gene files. |

### Usage

#### RankCompV3 package in julia

[RankCompV3.jl](https://github.com/pathint/RankCompV3.jl)

##### Use directly

```julia
using RankCompV3
#runing code, the default values shown below can be used accordingly.
reoa("expr.txt",
    "metadata.txt";
    expr_threshold = 0,
    min_profiles = 0,
    min_features = 0,
    pval_reo = 0.01,
     pval_deg = 1.00,
     padj_deg = 0.05,
    use_pseudobulk = 0,
    use_hk_genes = "yes"
    hk_file = "HK_genes_info.tsv",
    gene_name_type = "ENSEMBL",
    ref_gene_max = 3000,
    ref_gene_min = 100
    n_iter = 128,
    n_conv = 5,
    cell_drop_rate = 0,
    gene_drop_rate = 0,
    work_dir = "./",
    use_testdata = "no"
    )
```

##### Local run usage

```julia
#Load the RankCompV3 package in julia. [path] indicates the storage path.
include("[path]/RankCompV3/src/RankCompV3.jl")
using Main.RankCompV3
#runing code, the default values shown below can be used accordingly.
reoa("expr.txt",
    "metadata.txt";
    expr_threshold = 0,
    min_profiles = 0,
    min_features = 0,
    pval_reo = 0.01,
     pval_deg = 1.00,
     padj_deg = 0.05,
    use_pseudobulk = 0,
    use_hk_genes = "yes"
    hk_file = "HK_genes_info.tsv",
    gene_name_type = "ENSEMBL",
    ref_gene_max = 3000,
    ref_gene_min = 100
    n_iter = 128,
    n_conv = 5,
    cell_drop_rate = 0,
    gene_drop_rate = 0,
    work_dir = "./",
    use_testdata = "no"
    )
```

### Output File Description

#### Resulting file

##### 1、fn_expr_ctrl.tsv

This file is the expression profile of the ctrl group.

##### 2、fn_expr_treat.tsv

This file is the expression profile of the treat group.

##### 3、fn_expr_result.tsv

This file contains Name, pval, padj, n11, n12, n13, n21, n22, n23, n31, n32, n33, Δ1, Δ2, se, z1.

##### 4、fn_expr_expr_dist.pdf

This file displays a graph of Distribution of Expression Values.

##### 5、fn_expr_expr_heat.pdf

This file shows heat maps of expression values for the ctrl and treat groups.

##### 6、fn_expr_contigency_table.pdf

This file shows the distribution of parameters in 3 x 3 contingency tables.

##### 7、fn_expr_delta_value.pdf

This file shows the delta distribution.

##### 8、fn_expr_se.pdf

This file shows the distribution of Standard Error (SE).

##### 9、fn_expr_z1.pdf

This file shows the z1 distribution.

##### 10、fn_expr_p_value.pdf

This file shows the distribution of p and fdr values.

##### 11、fn_expr_degs_expr_dist.pdf

This file shows the distribution of expression values for DEGs.

##### 12、fn_expr_degs_expr_heat.pdf

This file shows a heat map of the expression values of DEGs in the ctrl and treat groups.

#### log file

### Demonstration

#### Application

##### RankCompV3 package in julia

[RankCompV3.jl](https://github.com/pathint/RankCompV3.jl)

###### Use directly

```julia
#For details about how to download the RankCompV3 package, see 4.
using RankCompV3
#The package comes with test data. Use the default parameters. If you need to modify the parameters, add them directly.
reoa(use_testdata="yes")
#Or local data. Use Default parameters. If you want to modify the parameters, add them directly.
reoa("/public/yanj/data/fn_expr.txt",
	"/public/yanj/data/fn_metadata.txt"
)
```

###### Local run usage

```julia
#For details about how to download the RankCompV3 package, see 4.
##[path] indicates the storage path.
include("[path]/RankCompV3/src/RankCompV3.jl")
using Main.RankCompV3
#The package comes with test data. Use the default parameters. If you need to modify the parameters, add them directly.
reoa(use_testdata="yes")
#Or local data. Use Default parameters. If you want to modify the parameters, add them directly.
reoa("/public/yanj/data/fn_expr.txt",
	"/public/yanj/data/fn_metadata.txt"
)
```

#### Input File

- Expression profile file  

 	[fn_expr.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_expr.txt)

- metadata file

 	[fn_metadata.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_metadata.txt)

- Housekeeping gene file (built-in, also supported for re-provisioning)

 	[HK_genes_info.tsv](https://github.com/yanjer/RankCompV3.jl/blob/master/hk_gene_file/HK_genes_info.tsv)

#### Resulting file

##### 1、fn_expr_ctrl.tsv

This file is the expression profile of the ctrl group.

```
[fn_expr_ctrl.tsv](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_ctrl.tsv)
```

##### 2、fn_expr_treat.tsv

This file is the expression profile of the treat group.

```
[fn_expr_treat.tsv](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_treat.tsv)
```

##### 3、fn_expr_result.tsv

This file contains Name, pval, padj, n11, n12, n13, n21, n22, n23, n31, n32, n33, Δ1, Δ2, se, z1.

```
[fn_expr_result.tsv](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_result.tsv)
```

##### 4、fn_expr_expr_dist.pdf

This file displays a graph of Distribution of Expression Values.

```
[fn_expr_expr_dist.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_expr_dist.pdf)
```

##### 5、fn_expr_expr_heat.pdf

This file shows heat maps of expression values for the ctrl and treat groups.

```
[fn_expr_expr_heat.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_expr_heat.pdf)
```

##### 6、fn_expr_contigency_table.pdf

This file shows the distribution of parameters in 3 x 3 contingency tables.

```
[fn_expr_contigency_table.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_contigency_table.pdf)
```

##### 7、fn_expr_delta_value.pdf

This file shows the delta distribution.

```
[fn_expr_delta_value.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_delta_value.pdf)
```

##### 8、fn_expr_se.pdf

This file shows the distribution of Standard Error (SE).

```
[fn_expr_se.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_se.pdf)
```

##### 9、fn_expr_z1.pdf

This file shows the z1 distribution.

```
[fn_expr_z1.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_z1.pdf)
```

##### 10、fn_expr_p_value.pdf

This file shows the distribution of p and fdr values.

```
[fn_expr_p_value.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_p_value.pdf)
```

##### 11、fn_expr_degs_expr_dist.pdf

This file shows the distribution of expression values for DEGs.

```
[fn_expr_degs_expr_dist.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_degs_expr_dist.pdf)
```

##### 12、fn_expr_degs_expr_heat.pdf

This file shows a heat map of the expression values of DEGs in the ctrl and treat groups.

```
[fn_expr_degs_expr_heat.pdf](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_degs_expr_heat.pdf)
```

#### log file

[RankCompV3-test-data-output.log](https://github.com/yanjer/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/RankCompV3_test_data_output.log)

### Suggestion and Issue Reporting

Any suggestion or issue reporting is welcome! You can contact **wang.xianlong@139.com** and **yanjer123@qq.com**. 
