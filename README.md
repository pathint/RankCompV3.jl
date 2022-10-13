## RankCompV3

​	RankCompV3 is a differentially expressed gene recognition algorithm based on relative expression order relation REO. The tool is developed based on the julia language, and the software is available for direct use. The details are described below. julia recommends using version 1.7 or later.

**RankCompV3 package in julia：https://github.com/Yanjj1/RankCompV3.jl.git**

**RankCompV3 software：https://github.com/Yanjj1/RankCompV3-software.git**

### Scope of application of RankCompV3

| Classification                              | Applicable scope                                             |
| ------------------------------------------- | ------------------------------------------------------------ |
| **Species**                                 | Unrestricted species                                         |
| **Species supported by housekeeping genes** | human/mouse                                                  |
| **Data type**                               | Microarray data, high-throughput sequencing data (RNA-seq), single-cell sequencing data (scRNA-seq), Cell expression by linear amplification and sequencing (CEL-seq) , methylation data and proteome data |
| **Data format**                             | Count/log2Count/FPKM/RPKM/TPM/Expression values              |
| **Genotype**                                | REFSEQ/SYMBOL/ENTREZID/ENSEMBL/UNIGENE/GENENAME              |

### Tool Dependency package

​	**Where julia related software packages have been added to the environment. The software packages related to the R language will be downloaded automatically, or you can download them manually.**

#### Package in julia:

```julia
    Distributed,
    SharedArrays,
    MultipleTesting (>= v0.5.1),
    DataFrames (>= v1.4.1),
    DelimitedFiles,
    CSV (>= v0.10.4),
    Statistics,
    RCall (>= v0.13.13),
    HypothesisTests (>= v0.10.10),
    Distributions (>= v0.25.75),
    Random,
    LinearAlgebra,
    ArgParse (>= v1.1.4)
```

#### Package in R:

```R
	ggplot2 (>= 3.3.6),
	pheatmap (>= 1.0.12)
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
Pkg.add(url="https://github.com/Yanjj1/RankCompV3.jl.git")
```

##### Local run usage

```shell
#configured in linux
#clone the RankCompV3 package from github to local
git clone https://github.com/Yanjj1/RankCompV3.jl.git
#Load the project dependency package
#[path] is the path from git clone https://github.com/Yanjj1/RankCompV3.git
julia --project=RankCompV3 [path]/RankCompV3/src/RankCompV3.jl
```

#### RankCompV3 software

[Yanjj1/RankCompV3-software](https://github.com/Yanjj1/RankCompV3-software)

```shell
#configured in linux
#clone RankCompV3 software from github to local
git clone https://github.com/Yanjj1/RankCompV3-software.git
#unzip
unzip RankCompV3-software.zip
```

### Parameters for details

| Parameter        | Parameter types              | Default value     | Parameters to describe                                       |
| ---------------- | ---------------------------- | ----------------- | ------------------------------------------------------------ |
| fn_expr          | AbstractString               | fn_expr.txt       | Gene expression profile file path. (required)                |
| fn_metadata      | AbstractStringAbstractString | fn_metadata.txt   | Grouping information file path. (required)                   |
| expr_threshold   | NumberNumber                 | 0                 | Gene expression threshold.                                   |
| pval_reo         | AbstractFloatAbstractFloat   | 0.01              | Stable threshold for p-value.                                |
| pval_sign_reo    | AbstractFloatAbstractFloat   | 1                 | Significant  reversal threshold for p-value.                 |
| padj_sign_reo    | AbstractFloatAbstractFloat   | 0.05              | Significant reversal threshold for FDR  value.               |
| hk_file          | AbstractString               | HK_genes_info.tsv | Housekeeper gene  file path.                                 |
| hk_name          | AbstractString               | ENSEMBL           | Column name of the  column where the housekeeping gene is located. |
| ref_gene_num     | Int                          | 3000              | The upper limit of  the number of housekeeping genes, if it is greater than this value,  ref_gene_num housekeeping genes are randomly selected from it. |
| use_housekeeping | AbstractString               | yes               | Do not use  housekeeping gene to set 1, use housekeeping gene to set 0. |
| species          | AbstractString               | human             | Data for species information, current housekeeping genes support only human and mouse species.Click [human](https://github.com/Yanjj1/RankCompV3/blob/master/RankCompV3/hk_gene_file/HK_genes_info.tsv) or [mouse](https://github.com/Yanjj1/RankCompV3/blob/master/RankCompV3/hk_gene_file/HK_genes_info_mouse.tsv) to see the housekeeping gene file. |
| cell_drop_rate   | Int                          | 0                 | At least how many genes were detected in each sample. (value non-zero). |
| gene_drop_rate   | Int                          | 0                 | At least how many cells each gene was  detected in (value non-zero). |
| work_dir         | AbstractString               | ./                | Working Directory.                                           |
| use_testdata     | AbstractString               | no                | Whether to use the default provided test data for analysis, yes or no. |

### Input File Description

| File                        | Description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| **Expression profile file** | Each row represents the gene, each column represents the sample and the expression matrix of the gene in the first column. |
| **Metadata file**           | First column sample name, second column group information.   |
| **housekeeping gene file**  | Housekeeping gene file with built-in support for several gene forms: Name, REFSEQ, SYMBOL, ENTREZID, ENSEMBL, UNIGENE and GENENAME. Also supports external incoming or not using housekeeping gene files. |

### Usage

#### RankCompV3 package in julia

[Yanjj1/RankCompV3.jl](https://github.com/Yanjj1/RankCompV3)

##### Use directly

```julia
using RankCompV3
#runing code
RankCompV3.reoa("expr.txt",
        "metadata.txt";
    expr_threshold = 0,
          pval_reo = 0.01,
     pval_sign_reo = 1.00,
     padj_sign_reo = 0.05,
           hk_file = "HK_genes_info.tsv",
           hk_name = "ENSEMBL",
      ref_gene_num = 3000,
    no_use_housekeeping = 0,
           species = "human",
    cell_drop_rate = 0,
    gene_drop_rate = 0,
    work_dir = "./",
    use_testdata = "no"
    )
```

##### Local run usage

```julia
#Load the RankCompV3 package in julia
include("[path]/RankCompV3/src/RankCompV3.jl")
using Main.RankCompV3
#runing code
reoa("expr.txt",
     "metadata.txt";
    expr_threshold = 0,
          pval_reo = 0.01,
     pval_sign_reo = 1.00,
     padj_sign_reo = 0.05,
           hk_file = "HK_genes_info.tsv",
           hk_name = "ENSEMBL",
      ref_gene_num = 3000,
    no_use_housekeeping = 0,
           species = "human",
    cell_drop_rate = 0,
    gene_drop_rate = 0,
    work_dir = "./",,
    use_testdata = "no"
    )
```

#### RankCompV3 software

[Yanjj1/RankCompV3-software](https://github.com/Yanjj1/RankCompV3-software)

```shell
#Used in linux
#See the help
RankCompV3-software/bin/RankCompV3 -h
#runing code
RankCompV3-software/bin/RankCompV3 --fn_expr "fn_expr.txt" --fn_metadata "fn_metadata.txt" --expr_threshold 0 --pval_reo 0.01 --pval_sign_reo 1.00 --padj_sign_reo 0.05 --hk_file "/home/yanj/jupyter_work/McCullagh/outcome_complete/HK_genes_info.tsv" --hk_name "ENSEMBL" --ref_gene_num 3000 --use_housekeeping "yes" --species "human" --cell_drop_rate 0 --gene_drop_rate 0 --work_dir ./ --use_testdata "no"
```

### Output File Description

#### Resulting file

##### With the housekeeping gene, the following seven files are generated. Otherwise, it is not generated.

- Three samples were randomly selected from each group to draw the expression distribution map of housekeeping gene in the samples.

- Gene expression profile files of labeled housekeeping genes, with rows representing genes and columns representing samples. The first column is the name of the gene, and the second column is whether it is a housekeeping gene.

##### Each iteration produces six files.

- The differential gene result file obtained in this iteration. Among them, there are 5 columns, each column is gene name, delta, sd_delta, p, FDR.

- The distribution maps of delta, sd_delta, Pval and Padj of the whole genes.

- The gene expression profile file of the marker housekeeping gene in this iteration, with rows representing genes and columns representing samples. The first column is the name of the gene, and the second column is whether it is a housekeeping gene.

##### Final output file after differential gene stabilization.

- The final calculation results (gene name, delta, sd_delta, p, FDR) without threshold screening.

- Differential gene expression profile.

- Differential gene expression heat map file.

#### log file

### Application example

#### code

#### RankCompV3 package in julia

[Yanjj1/RankCompV3.jl](https://github.com/Yanjj1/RankCompV3)

##### Use directly

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

##### Local run usage

```julia
#For details about how to download the RankCompV3 package, see 4.
using Main.RankCompV3
#The package comes with test data. Use the default parameters. If you need to modify the parameters, add them directly.
reoa(use_testdata="yes")
#Or local data. Use Default parameters. If you want to modify the parameters, add them directly.
reoa("/public/yanj/data/fn_expr.txt",
	"/public/yanj/data/fn_metadata.txt"
)
```

#### RankCompV3 software

[Yanjj1/RankCompV3-software](https://github.com/Yanjj1/RankCompV3-software)

```shell
#See the help
RankCompV3/bin/RankCompV3 -h
#The package comes with test data. Use the default parameters. If you need to modify the parameters, add them directly.
RankCompV3/bin/RankCompV3 --use_testdata "yes"
#Or local file. Use Default parameters. If you want to modify the parameters, add them directly.
RankCompV3/bin/RankCompV3 --fn_expr "/public/yanj/data/fn_expr.txt" --fn_metadata "/public/yanj/data/fn_metadata.txt"
```

#### Input File

- Expression profile file  

​			[fn_expr.txt](https://github.com/Yanjj1/RankCompV3.jl/blob/master/RankCompV3/test/fn_expr.txt)

- metadata file

​			[fn_metadata.txt](https://github.com/Yanjj1/RankCompV3.jl/blob/master/RankCompV3/test/fn_metadata.txt)

- Housekeeping gene file (built-in, also supported for re-provisioning)

​			[HK_genes_info.tsv](https://github.com/Yanjj1/RankCompV3.jl/blob/master/RankCompV3/hk_gene_file/HK_genes_info.tsv)

#### Resulting file

##### With the housekeeping gene, the following seven files are generated. Otherwise, it is not generated.

- Three samples were randomly selected from each group to draw the expression distribution map of housekeeping gene in the samples.

​			[fn_expr_hk_nonhk_gene_sample1.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_hk_nonhk_gene_sample1.pdf)

​			[fn_expr_hk_nonhk_gene_sample2.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV	3-test-data-output/fn_expr_hk_nonhk_gene_sample2.pdf)

​			[fn_expr_hk_nonhk_gene_sample4.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_hk_nonhk_gene_sample4.pdf)

​			[fn_expr_hk_nonhk_gene_sample5.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_hk_nonhk_gene_sample5.pdf)

​			[fn_expr_hk_nonhk_gene_sample6.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_hk_nonhk_gene_sample6.pdf)

​			[fn_expr_hk_nonhk_gene_sample8.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_hk_nonhk_gene_sample8.pdf)

- Gene expression profile files of labeled housekeeping genes, with rows representing genes and columns representing samples. The first column is the name of the gene, and the second column is whether it is a housekeeping gene.

​			[fn_expr_hk_nonhk_gene.tsv](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_hk_nonhk_gene.tsv)

##### Each iteration produces six files. (Only the results of the 0th iteration are shown below)

- The differential gene result file obtained in this iteration. Among them, there are 5 columns, each column is gene name, delta, sd_delta, p, FDR.

​			[fn_expr_iteration_0_result.tsv](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_iteration_0_result.tsv)

- The distribution maps of delta, sd_delta, Pval and Padj of the whole genes.

​		delta: 

​			[fn_expr_0_delta_graph.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_0_delta_graph.pdf)

​		sd_delta: 

​			[fn_expr_0_sd_delta_graph.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_0_sd_delta_graph.pdf)

​		*p*-value: 

​			[fn_expr_0_Pval_graph.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_0_Pval_graph.pdf)

​		FDR: 

​			[fn_expr_0_Padj_graph.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_0_Padj_graph.pdf)

- The gene expression profile file of the marker housekeeping gene in this iteration, with rows representing genes and columns representing samples. The first column is the name of the gene, and the second column is whether it is a housekeeping gene.

​			[fn_expr_hk_nonhk_gene_0.tsv](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_hk_nonhk_gene_0.tsv)

##### Final output file after differential gene stabilization.

- The final calculation results (gene name, delta, sd_delta, p, FDR) without threshold screening.

​			[fn_expr_deg_exp.tsv](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_deg_exp.tsv)

- Differential gene expression profile.

​			[fn_expr_result.tsv](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_result.tsv)

- Differential gene expression heat map file.

​			[fn_expr_deg_exp_graph.pdf](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/fn_expr_deg_exp_graph.pdf)

#### log file

​	[RankCompV3-test-data-output.log](https://github.com/Yanjj1/RankCompV3-test-data-output/blob/master/RankCompV3-test-data-output/RankCompV3-test-data-output.log)

### Suggestion and Issue Reporting

​	Any suggestion or issue reporting is welcome! You can contact yanjer123@qq.com. 
