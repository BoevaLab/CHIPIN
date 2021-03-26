



# CHIPIN
ChIP-seq Intersample Normalization

Copywrite Lélia Polit, 2020

This program is free software; you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
 

## INTRODUCTION

CHIPIN is an R package that provides a normalization procedure to compare signals between ChIP-seq samples in different conditions; experiments should be performed using the same antibody. CHIPIN output includes normalized density profiles and several statistics describing the characteristics of the normalization procedure. The normalization is based on the assumption that no differences in ChIP-seq signal should be observed in the regulatory regions of genes whose expression does not change across samples/conditions. Moreover, CHIPIN provides a possibility to qualify the antibody used in the ChIP-seq experiments by plotting ChIP-seq signal density around gene transcription start sites for highly, medium and lowly expressed genes.


## REQUIREMENTS

1. The R libraries listed below are required. They should be *automatically* installed when you run the intallation of CHIPIN via devtools:


```
install.packages(devtools)
library(devtools)
devtools::install_github("BoevaLab/CHIPIN", build_vignettes = TRUE)
```

You can also install them manually:
CRAN libraries:
```
install.packages(c("stringr", "ggplot2", "matrixStats", "viridis", "tidyr", "plotrix", "rlist", "RColorBrewer", "tiger", "pracma", "dplyr", "gridExtra", "rminer"))
```
Bioconductor libraries:
```
BiocManager::install(c("rtracklayer", "preprocessCore", "GenomicFeatures", "biomaRt", "ggpubr"))
```



2. Deeptools installation is required.
Please see additional Python libraries required for Deeptools.
https://deeptools.readthedocs.io/en/develop/content/installation.html


We recommend to create an `anaconda` environment in which Deeptools is installed. In this case you can launch RStudio from this environment and all the path variables are set correctly.

```
# create conda environment
conda create --name CHIPIN
conda activate CHIPIN

# install Deeptools
conda install -c bioconda deeptools

# launch RStudio from the terminal
open -na Rstudio
```

## INSTALLATION

Use devtools for installation of CHIPIN from GitHub:

```
install.packages(devtools)
library(devtools)
devtools::install_github("BoevaLab/CHIPIN", build_vignettes = TRUE)
```

Alternatively, after downloading the .zip archive, one should un-zip it and run: 
R CMD INSTALL pathTo/CHIPIN-master

## FUNCTIONS

### plot_expression function: 
CHIPIN offers the possibility to profile ChIP-seq intensity around TSS as a function of gene expression level using the function plot_expression. The signal is visualized for three groups of genes obtained with k-means clustering: highly-expressed, medium-expressed and lowly-expressed genes. The results of such visualization are important to verify the efficiency of the antibody used. To use this function, several parameters are mandatory:

* **TPM**: path to a gene expression file (TPM values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". If you provide the "TPM" parameter, do not use the "raw_read_count" parameter or the "RPKM" parameter. Default: NULL.
* **RPKM**: path to a gene expression file (RPKM values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". If you provide the "RPKM" parameter, do not use the "raw_read_count" parameter or the "TPM" parameter. Default: NULL.
* **raw_read_count**: path to a gene expression file (raw read count values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". Raw read count values will be transformed into RPKM values using information on exon lengths. If you provide the "raw_read_count" parameter, do not use the "RPKM" parameter or the "TPM" parameter. Default: NULL.
* **path_to_bw**: a vector containing paths to .bigWig files of the samples/conditions of interest. ! Mandatory parameter.
* **output_dir**: path to the output directory where one wants to store the ouput files. This directory should be created before running the function. Default: "." 
* **organism**: is one of “hg19”, “hg38”, “mm9” or “mm10”. ! Mandatory parameter with no default value. 

One optional parameter can be set:

* **histone_mark** name of the histone mark of interest; used to plot legends. Default:"ChIP-seq signal"

```
#### Usage:

##### using TPM values:
plot_expression(TPM, RPKM=NULL, raw_read_count=NULL, path_to_bw, output_dir=".", organism, histone_mark="ChIP-seq signal")

##### using RPKM values:
plot_expression(TPM=NULL, RPKM, raw_read_count=NULL, path_to_bw, output_dir=".", organism, histone_mark="ChIP-seq signal")

##### using raw read count values:
plot_expression(TPM=NULL, RPKM=NULL, raw_read_count, path_to_bw, output_dir=".", organism, histone_mark="ChIP-seq signal")
```

This function can be also launched directly in the CHIPIN_normalize function by setting the parameter expression_plot to TRUE.



### CHIPIN_normalize function
This is the main function of the package; it identifies genes that do not change their expression across the conditions ("constant_genes") and performs the normalization. All the parameters of this function are explained in the following sections. There are three steps performed by the CHIPIN_normalize function:

* Determine constant_genes
* Perform normalization
* Compute statistics

Common parameters for the three steps of CHIPIN_normalize function:

* **sample_name**: will be used as a prefix for the different outputs. Default:"sample"
* **output_dir**: path to the output directory where one wants to store the ouput files. This directory should be created before running the function. Default: "."
* **organism**: is one of “hg19”, “hg38”, “mm9” or “mm10”.

To determine “constant_genes”, the mandatory parameters are:

* **TPM**: path to a gene expression file (TPM values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". TPM values will be transformed into "raw_read_count" values using information on exon lengths; then, "raw_read_count" values will be used to determine genes whose expression does not change across all the conditions ("constant_genes"). If you provide the "TPM" parameter, do not use the "raw_read_count" parameter or the "RPKKM" parameter. If both "RPKM", "raw_read_count" and "TPM" parameters are set to NULL, and "path_to_file_with_constant_genes" is NULL too, then all genes will be used for the normalization; "expression_plot" (see below) will be set to FALSE. Default: NULL
* **RPKM**: path to a gene expression file (RPKM values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". RPKM values will be transformed into "raw_read_count" values using information on exon lengths; then, "raw_read_count" values will be used to determine genes whose expression does not change across all the conditions ("constant_genes"). If you provide the "RPKM" parameter, do not use the "raw_read_count" parameter or the "TPM" parameter. If both "RPKM", "raw_read_count" and "TPM" parameters are set to NULL, and "path_to_file_with_constant_genes" is NULL too, then all genes will be used for the normalization; "expression_plot" (see below) will be set to FALSE. Default: NULL
* **raw_read_count**: path to a gene expression file (raw read count values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". If you provide the "raw_read_count" parameter, do not use the "RPKM" parameter or the "TPM" parameter. If both "RPKM", "raw_read_count" and "TPM" parameters are set to NULL, and "path_to_file_with_constant_genes" is NULL too, then all genes will be used for the normalization; "expression_plot" (see below) will be set to FALSE. Default: NULL

Optional paramaters:
* **percentage**: a value between 0 and 1 describing the percentage of the total number of genes that one wants to be defined as "constant_genes". Default: 0.1
* **path_to_file_with_constant_genes**: path to a .bed file with genes that do not change their expression across the conditions ("constant_genes"). If left emtpy (NULL), the list of constant genes will be determined automatically using either "RPKM" or "raw_read_count" values. Default:NULL

To perform normalization, the mandatory parameter is:

* **path_to_bw**: a vector containing paths to .bigWig files of the samples/conditions of interest. ! Mandatory parameter with no default value

Optional parameters:
* **type_norm**: type of normalization to perform: 'linear' or 'quantile'. Default: 'linear'
* **beforeRegionStartLength** (default: 4000), **afterRegionStartLength** (default: 4000), **regionBodyLength** (default: 40000), **binSize** (default: 10): parameters of the “computematrix” function of deeptools. They correspond to distance upstream of the reference-point selected, distance downstream of the reference-point selected, distance in bases to which all regions will be fit, and length, in bases, of the non-overlapping bins for averaging the score over the regions length, respectively. See https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details
* **expression_plot**: boolean parameter, use "expression_plot=TRUE"" to call function “plot_expression” to plot the density signal around gene TSS. Default: FALSE
* **compute_stat**: boolean parameter, use "compute_stat=TRUE" to compute statisctics characterizing the normalization process. This statistic will be written in the "output_StatsFile.txt" file located in the output_folder and will show how much the normalization reduced the difference between the samples/conditions. Default:FALSE
* **nGroup**: number of gene groups for quantile normalisation. Default: 20
* **histone_mark**: name of the histone mark of interest; used to plot legends. Default:"ChIP-seq signal"
* **nThreads**: number of processors to use. Default: 1

```
#### Usage:
CHIPIN_normalize(path_to_bw, ...)

##### Using TPM values:
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM, RPKM=NULL, raw_read_count=NULL, path_to_file_with_constant_genes=NULL, sample_name, output_dir=".", organism, beforeRegionStartLength=4000, afterRegionStartLength=4000, regionBodyLength=40000, binSize=10, expression_plot=FALSE, compute_stat=FALSE, percentage=0.1, nGroup=20, histone_mark="ChIP-seq signal", nThreads=1)

##### Using RPKM values:
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM=NULL, RPKM, raw_read_count=NULL, path_to_file_with_constant_genes=NULL, sample_name, output_dir=".", organism, beforeRegionStartLength=4000, afterRegionStartLength=4000, regionBodyLength=40000, binSize=10, expression_plot=FALSE, compute_stat=FALSE, percentage=0.1, nGroup=20, histone_mark="ChIP-seq signal", nThreads=1)


##### Using raw read count values:
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM=NULL, RPKM=NULL, raw_read_count, path_to_file_with_constant_genes=NULL, sample_name, output_dir=".", organism, beforeRegionStartLength=4000, afterRegionStartLength=4000, regionBodyLength=40000, binSize=10, expression_plot=FALSE, compute_stat=FALSE, percentage=0.1, nGroup=20, histone_mark="ChIP-seq signal", nThreads=1)


##### Using constant genes provided by the user:
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM=NULL, RPKM=NULL, raw_read_count=NULL, path_to_file_with_constant_genes, sample_name, output_dir=".", organism, beforeRegionStartLength=4000, afterRegionStartLength=4000, regionBodyLength=40000, binSize=10, expression_plot=FALSE, compute_stat=FALSE, percentage=0.1, nGroup=20, histone_mark="ChIP-seq signal", nThreads=1)


##### Using all genes (not recommended):
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM=NULL, RPKM=NULL, raw_read_count=NULL, path_to_file_with_constant_genes=NULL, sample_name, output_dir=".", organism, beforeRegionStartLength=4000, afterRegionStartLength=4000, regionBodyLength=40000, binSize=10, expression_plot=FALSE, compute_stat=FALSE, percentage=0.1, nGroup=20, histone_mark="ChIP-seq signal", nThreads=1)


```


## EXAMPLES: 

```
library(CHIPIN)

#initialize parameters:
pathToRPKMfile = system.file("extdata", "FPKM_values_CLBBER_CLBMA_SJNB12.txt", package = "CHIPIN")
pathToFiles = system.file("extdata", c("CLBBER.K27ac.rep3.bw","SJNB12.K27ac.rep3.bw","CLBMA.K27ac.rep3.bw"), package = "CHIPIN")
outputFolder = "." #create the corresponding output folder if it does not exists
histoneMarkName = "H3K27Ac"
sampleName = "neuroblastoma"


#normalize the data without plotting the distribution around gene TSS (quantile normalization, expression_plot=FALSE):
CHIPIN_normalize(path_to_bw=pathToFiles, type_norm="quantile", RPKM=pathToRPKMfile, sample_name=sampleName, output_dir=outputFolder, organism="hg19", compute_stat=TRUE, percentage=0.1, nGroup=20, histone_mark=histoneMarkName)


#normalize the data and plot the distribution around gene TSS (linear normalization, expression_plot=TRUE):
CHIPIN_normalize(path_to_bw=pathToFiles, type_norm="linear", RPKM=pathToRPKMfile, sample_name=sampleName, output_dir=outputFolder, organism="hg19", expression_plot=TRUE, compute_stat=TRUE, histone_mark=histoneMarkName)


```

Examples for running the two functions are also provided in the vignette file included in the package: file "CHIPIN-vignette.html". To open the HTML vignette, you can type:

```
openVignette("CHIPIN")
```


