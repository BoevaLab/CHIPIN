#Introduction

CHIPIN is a R package that aims to provide a normalization procedure to compare ChIP-seq signals between different conditions/samples when spike-in information is not available. The normalization procedure is based on the assumption that on average, no difference in ChIP-seq signal should be observed in the regulatory regions of genes whose expression levels are constant across samples/conditions. Using CHIPIN, the user obtains bigWig files that can be further used to compare intensity of histone modifications or transcription factor binding at specific regions of interest between the different conditions/samples.


```{r setup}
library(CHIPIN)
```


#ChIPIN_normalize function

This section will explain how the normalization procedure is done. To perform the normalization, the main function: ChIPIN_normalize should be used. All the parameters of this function will be explained in the following sections. At the end of this section, there will be an example of how to use this function with the test data included in the package. There is three parts in the ChIPIN_normalize function:

+ Determine ConstantGenes
+ Perform normalization
+ Compute statistics

##Determine ConstantGenes

If the user provides a list of genes in bed format through a bed file in the **pathToFileGenesNotMoving** parameter, the package will use this set of genes for the normalization. 

On the contrary, the parameter **pathToFileGenesNotMoving** can be set to NULL, and ChIPIN will determine the set of genes whose expression levels are constant across samples/conditions. To this end, one should provide RPKM or raw_read_count data. Given the gene expression data (RNA-seq or microarray), the mean and the standard deviation of the Count Per Million (CPM) values for each gene across samples/conditions will be determined. Ten percent (this percentage can be modified using the parameter percentage) of genes showing the smallest standard deviation is extracted as “ConstantGenes”. The output is a standard bed file that will be stored in the directory the user provides in the “output_dir” paramater. 

However, if the user can not provide either RPKM/raw_read_count data nor list of "ConstantGenes", all genes will further be used in the normalization process eventhough  we do not recommend this way as it may result in improper correction.



##Perform Normalization

The normalization process starts by building the matrix containing the original binding intensities of ChIP-seq signal across “ConstantGenes” (regions +/- x bp surrounding the gene body, default: x=4kb) using the function “computematrix” included in the deeptools package (Ramirez et al, 2016, 2).

These matrices obtained for each sample/condition are used to infer the normalization parameters. Then two different types of normalization can be performed:

+ **Quantile normalization**
+ **Linear regression with non-zero intercept**

![](CHIPIN/vignettes/figS2-page-001.jpg){width=85%}


###Quantile normalization

In order to perform quantile normalization, the parameter **typeNorm** should be set to **“quantile”**. The output matrix of the function “computematrix” is sorted by rows according to the total signal value of each gene across bins. Given the ranked matrix, we built k groups (k=20 by default) corresponding to k different ChIP-seq signal intensities. The quantile normalization is done on the mean values of these k groups. The main steps of the quantile normalization are explained on the following figure.


![](fig2-page-001.jpeg){width=85%}

###Linear normalization with non-zero intercept

For the linear normalization with non-zero intercept, the parameter **typeNorm** should be set to **“linear”**. The linear regression with non-zero intercept is performed on the average signal intensity values per bin i for a given sample versus a reference sample. The reference sample is choosen as the sample with the median total signal value among all samples. For each sample, the parameters α and β that minimize the sum of square errors are determined by the linear regression with non-zero intercept. The signal of the bigWig file for the current sample is then modified using these parameters α and β.

##Compute statistics

CHIPIN computes statistics illustrating success of the normalization process: the relative difference between average signal curves before and after the normalization. For each sample, CHIPIN computes the area under the average signal curves before the normalization around TSS of “ConstantGenes”. The sample with the highest value is selected as “reference”. Then the value of each other sample is expressed as a percentage of the “reference”. Given these values, CHIPIN computes the relative difference between all samples before normalization. After the normalization, the same process is done. If the normalization process worked well the relative difference between the percentage of all samples should decrease after the normalization.


##Parameters

Common parameters for the three parts of ChIPIN_normalize function:

+ **sample_name** will be used as a prefix for the different parameters.
+ **output_dir** is the path to the directory where one wants to store all output of the ChIPIN_normalize function. This directory should be created before running the function.
+ **organism** is one of “hg19”, “hg38”, “mm9” or “mm10”.

To determine "ConstantGenes", the mandatory parameters are:

+ **RPKM** is the path to RPKM values, each column of the data should correspond to one sample/condition. If one provides raw_read_count data, put NULL for raw_read_count parameter.
+ **raw_read_count** is the path to raw read count values, each column of the data should correspond to one sample/condition. If one provides raw_read_count data, put NULL for RPKM parameter.

Other paramaters can be modified :

+ **percentage** (default: 10) corresponds to the percentage of all genes that will be selected as ConstantGenes.
+ **pathToFileGenesNotMoving** (default: NULL) is the path to a bed file containing "ConstantGenes".

To perform normalization, mandatory parameters are: 

+ **fileWithPaths** is a textfile with lines corresponding to the paths to the bigWig files.

Additional parameters can be modified:

+ **beforeRegionStartLength** (default: 4000), **afterRegionStartLength** (default: 4000), **regionBodyLength** (default: 40000), **binSize** (default: 10) are “computematrix” function’s parameters.
+ **expression_plot** (default: FALSE) is a boolean parameter, if equal to TRUE the function plot_expression, that wil be further explained in this vignette, will be run.
+ **compute_stat** (default: FALSE) is a boolean parameter, if equal to TRUE, CHIPIN will compute statistics illustrating sucess of the normalization process (*i. e.* the relative difference between average signal curves before and after the normalization).
+ **typeNorm** (default: “quantile”) can be “quantile” or “linear”.
+ **nGroup** (default: 20) is a parameter dedicated to quantile normalization.
+ **histoneMark** (default: "ChIP-seq signal") is used in the legend of plots.


```
ChIPIN_normalize(RPKM, raw_read_count, pathToFileGenesNotMoving=NULL, sample_name, output_dir, organism, fileWithPaths, beforeRegionStartLength=4000, afterRegionStartLength=4000, regionBodyLength=40000, binSize=10, expression_plot=FALSE, compute_stat=FALSE, percentage=10, typeNorm="quantile", nGroup=20, histoneMark="ChIP-seq signal")

```


##Example

Using the test data included in the package, to perform quantile normalization, the parameters are:

+ **RPKM**: system.file("data", "FPKM_values_CLBBER_CLBMA_SJNB12.txt", package = "CHIPIN") 
+ **raw_read_count**: NULL
+ **pathToFileGenesNotMoving**: NULL
+ **sample_name**: “neuroblastoma”
+ **output_dir**: "/.../test_ChIPIN/" **(this directory should be previously created.)**
+ **organism**: “hg19”
+ **fileWithPaths**: paste(output_dir, "pathToFiles2.txt", sep="")
+ **beforeRegionStartLength**: 4000
+ **afterRegionStartLength**: 4000
+ **regionBodyLength**: 40000
+ **binSize**: 10
+ **expression_plot**: FALSE
+ **compute_stat**: TRUE
+ **percentage**: 10
+ **typeNorm**: “quantile”
+ **nGroup**: 20
+ **histoneMark** : "H3K27Ac"

Before running the function **ChIPIN_normalize**, there are two things one should do: 

+ First, create the **output_dir** directory. Then the path to **output_dir** should be provided in the following format: "/path/To/Directory**/**". Please do not forget the "/" at the end of the path.  
+ Second, create the file containing the paths to the bigwig files of the test data. In case of the test data, run this command: 

```
write(c(system.file("data", "CLBBER.K27ac.rep3.bw", package = "CHIPIN"), system.file("data", "SJNB12.K27ac.rep3.bw", package = "CHIPIN"), system.file("data", "CLBMA.K27ac.rep3.bw", package = "CHIPIN")), file=paste(output_dir, "pathToFiles2.txt", sep=""))

```

The file should look like this:

![](paths.png){width=85%}

Then to execute the function ChIPIN_normalize, run this command:

```
ChIPIN_normalize(system.file("data", "FPKM_values_CLBBER_CLBMA_SJNB12.txt", package = "CHIPIN"), NULL, NULL, "neuroblastoma", output_dir, "hg19", paste(output_dir, "pathToFiles2.txt", sep=""), 4000, 4000, 40000, 10, FALSE, TRUE, 10, "quantile", 20, "H3K27Ac")

```

##Output files

There is 9 different type of output files: 

+ **sample_nameGenesNotMoving.bed** is the bed file of "Constant Genes".

![](bedGenesNotmoving.png){width=50%}

+ **sample_nameCPMmeanVSsd.pdf** figure representing log2(standardDeviation(CPM)+1) function of log2(mean(CPM)+1) showing in red **"Constant genes"**.

![](neuroblastomaCPMmeanVSsd.png){width=50%}
                                                                                                                                               
+ **pathRenorm.txt** is a text file containing the path to the bigWig files normalized.
+ **normalized bigWig files** in case of linear regression the suffixe is "renormReg.bw", in case of quantile normalization the prefix is "QuantileNorm_".

![](igv_screen.png){width=100%}


+ **.mat.gz files** are the output of computeMatrix function from deeptools.
+ **StatsAfter.txt and StatsBefore.txt** are text files that are used to compute the statistics. 
+ **After_Normalization.pdf and Before_Normalization.pdf** are density profiles around TSS of "Constant Genes" that show the impact of normalization.

![](Before_Normalisation-1.png){width=50%}

![](After_NormalisationQuantile.png){width=50%}

![](After_NormalisationLinear.png){width=50%}



#plot_expression function

CHIPIN offers the possibility to profile ChIP-seq intensity around TSS as a function of gene expression level using the function plot_expression. The signal is visualized for three groups of genes obtained with k-means clustering: highly-expressed, medium-expressed and lowly-expressed genes. The results of such visualization are important to verify the efficiency of the antibody used. To use this function, several parameters are mandatory:

+ **RPKM** is the path to RPKM values, each column of the data should correspond to one sample/condition. If one provides raw_read_count data, put NULL for raw_read_count parameter.
+ **raw_read_count** is the path to raw read count values, each column of the data should correspond to one sample/condition. If one provides raw_read_count data, put NULL for RPKM parameter.
+ **fileWithPaths** is a textfile with lines corresponding to the paths to the bigWig files.
+ **output_dir** is the path to the directory where one wants to store all output of the ChIPIN_normalize function. This directory should be created before running the function.
+ **organism** is one of “hg19”, “hg38”, “mm9” or “mm10”.

One additional parameter can be modified:

+ **histoneMark** (default: "ChIP-seq signal") will be used for the legend of the figures.

```
plot_expression(RPKM, raw_read_count, fileWithPaths, output_dir, organism, histoneMark="ChIP-seq signal")
```


This function can be also launch directly in the **ChIPIN_normalize** function by setting the parameter **expression_plot** to TRUE.

##Example

Using the test data included in the package, the way to use the function plot_expression is:

+ **RPKM**: system.file("data", "FPKM_values_CLBBER_CLBMA_SJNB12.txt", package = "CHIPIN")
+ **raw_read_count**: NULL
+ **fileWithPaths**: system.file("data", "pathToFiles2.txt", package = "CHIPIN")
+ **output_dir**: "/.../test_ChIPIN/" (this directory should be previously created.)
+ **organism**: “hg19”
+ **histoneMark**: “H3K27Ac”


```
plot_expression(system.file("data", "FPKM_values_CLBBER_CLBMA_SJNB12.txt", package = "CHIPIN"), NULL, system.file("data", "pathToFiles2.txt", package = "CHIPIN"), output_dir, “hg19”, “H3K27Ac”)
```


![](expression-1.png){width=90%}

