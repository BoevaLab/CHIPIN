\name{CHIPIN_normalize}
\alias{CHIPIN_normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CHIPIN_normalize
%%  ~~function to do ... ~~
}
\description{
This is the main function of the package that should be used to find the “constant genes” and perform the normalization process.
}
\usage{
CHIPIN_normalize(...)

# Using TPM values:
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM, RPKM=NULL, raw_read_count=NULL, output_dir=".", organism, histone_mark="ChIP-seq signal")

#using RPKM values:
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM=NULL, RPKM, raw_read_count=NULL, output_dir=".", organism, histone_mark="ChIP-seq signal")

#using raw read count values:
CHIPIN_normalize(path_to_bw, type_norm="linear", TPM=NULL, RPKM=NULL, raw_read_count, output_dir=".", organism, histone_mark="ChIP-seq signal")
}

\arguments{
  \item{path_to_bw}{a vector containing paths to .bigWig files of the samples/conditions of interest. ! Mandatory parameter with no default value}
  \item{type_norm}{type of normalization to perform: 'linear' or 'quantile'. Default: 'linear'}
  \item{TPM}{path to a gene expression file (TPM values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". If you provide the "TPM" parameter, do not use the "raw_read_count" parameter or the "RPKM" parameter. If all gene expression parameters are set to NULL, and "path_to_file_with_constant_genes" is NULL too, then all genes will be used for the normalization; "expression_plot" (see below) will be set to FALSE. Default: NULL.}
  \item{RPKM}{path to a gene expression file (RPKM values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". RPKM values will be transformed into "raw_read_count" values using information on exon lengths; then, "raw_read_count" values will be used to determine genes whose expression does not change across all the conditions ("constant_genes"). If you provide the "RPKM" parameter, do not use the "raw_read_count" parameter or the "TPM" parameter. If all gene expression parameters are set to NULL, and "path_to_file_with_constant_genes" is NULL too, then all genes will be used for the normalization; "expression_plot" (see below) will be set to FALSE. Default: NULL}
  \item{raw_read_count}{path to a gene expression file (raw read count values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". If you provide the "raw_read_count" parameter, do not use the "RPKM" parameter or the "TPM" parameter. If all gene expression parameters are set to NULL, and "path_to_file_with_constant_genes" is NULL too, then all genes will be used for the normalization; "expression_plot" (see below) will be set to FALSE. Default: NULL}
  \item{path_to_file_with_constant_genes}{path to a .bed file with genes that do not change their expression across the conditions ("constant_genes"). If left emtpy (NULL), the list of constant genes will be determined automatically using either "RPKM" or "raw_read_count" values. Default:NULL}
  \item{sample_name}{sample name. Default: "sample"}
  \item{output_dir}{path to the output directory where one wants to store the ouput files. This directory should be created before running the function. Default: "."}
  \item{organism}{reference genome: "mm10", "mm9", "hg38" or "hg19". ! Mandatory parameter with no default value}
  \item{beforeRegionStartLength}{distance upstream of the reference-point selected - computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details. Default: 4000}
  \item{afterRegionStartLength}{distance downstream of the reference-point selected - computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details. Default: 4000}
  \item{regionBodyLength}{distance in bases to which all regions will be fit - computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details. Default: 40000}
  \item{binSize}{length, in bases, of the non-overlapping bins for averaging the score over the regions length - computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details. Default: 10}
  \item{expression_plot}{boolean parameter, use "expression_plot=TRUE"" to call function “plot_expression” to plot the density signal around gene TSS. Default: FALSE}
  \item{compute_stat}{boolean parameter, use "compute_stat=TRUE" to compute statisctics characterizing the normalization process. This statistic will be written in the "output_StatsFile.txt" file located in the output_folder and will show how much the normalization reduced the difference between the samples/conditions. Default: FALSE}
  \item{percentage}{a value between 0 and 1 describing the percentage of the total number of genes that one wants to be defined as "constant_genes". Default: 0.1}
  \item{nGroup}{number of gene groups for quantile normalisation. Default: 20}
  \item{histone_mark}{name of the histone mark of interest; used to plot legends. Default:"ChIP-seq signal"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lélia Polit, BoevaLab, "Computational Epigenetics of Cancer", Inserm, CNRS, Cochin Institute, Paris, France
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#initialize parameters:
pathToRPKMfile = system.file("extdata", "FPKM_values_CLBBER_CLBMA_SJNB12.txt", package = "CHIPIN")
pathToFiles = system.file("extdata", c("CLBBER.K27ac.rep3.bw","SJNB12.K27ac.rep3.bw","CLBMA.K27ac.rep3.bw"), package = "CHIPIN")
outputFolder = "." #change it if needed; create the corresponding output folder if it does not exists
histoneMarkName = "H3K27Ac"
sampleName = "neuroblastoma"


#normalize the data without plotting the distribution around gene TSS (quantile normalization, expression_plot=FALSE):
CHIPIN_normalize(path_to_bw=pathToFiles, type_norm="quantile", RPKM=pathToRPKMfile, sample_name=sampleName, output_dir=outputFolder, organism="hg19", compute_stat=TRUE, percentage=0.1, nGroup=20, histone_mark=histoneMarkName)


#normalize the data and plot the distribution around gene TSS (linear normalization, expression_plot=TRUE):
CHIPIN_normalize(path_to_bw=pathToFiles, type_norm="linear", RPKM=pathToRPKMfile, sample_name=sampleName, output_dir=outputFolder, organism="hg19", expression_plot=TRUE, compute_stat=TRUE, histone_mark=histoneMarkName)

}
