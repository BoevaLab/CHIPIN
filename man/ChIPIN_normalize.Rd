\name{ChIPIN_normalize}
\alias{ChIPIN_normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ChIPIN_normalize
%%  ~~function to do ... ~~
}
\description{
This is the main function of the package that should be used to find the “ConstantGenes” and to perform the normalization process.
}
\usage{

}

\arguments{
  \item{RPKM}{dataframe with two columns: "geneName", "RPKM". Will be transformed into raw_read_count values using exons length information, then raw_read_count values will be used to determine ConstantGenes. If you provide RPKM parameter, put NULL for raw_read_count parameter. default:NULL}
  \item{raw_read_count}{dataframe with two columns: "geneName", "raw_read_count". If you provide raw_read_count parameter, put NULL for RPKM parameter. default:NULL}
  \item{pathToFileGenesNotMoving}{bed file of ConstantGenes default:NULL}
  \item{sample_name}{sample name. default:"sample"}
  \item{output_dir}{path to output directory where one wants to store the ouput files. ! mandatory parameter}
  \item{organism}{can be mm10, mm9, hg38 or hg19. ! mandatory parameter}
  \item{fileWithPaths}{text file containing path to .bigWig files of the samples/conditions of interest. ! mandatory parameter}
  \item{beforeRegionStartLength}{computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details default: 4000}
  \item{afterRegionStartLength}{computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details default: 4000}
  \item{regionBodyLength}{computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details default: 40000}
  \item{binSize}{computeMatrix function parameter see https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html for more details default: 10}
  \item{expression_plot}{boolean parameter, TRUE means that the function “plot_expression” will be performed. default: FALSE}
  \item{compute_stat}{boolean parameter, TRUE means that a statistic will be computed after the normalization process. This statistic will be in the "output_StatsFile.txt" file in the output_folder and will show how much the normalization reduced the difference between the samples/conditions. default:FALSE}
  \item{percentage}{percentage of the total number of genes one wants as "ConstantGenes". default: 10}
  \item{typeNorm}{'linear' or 'quantile' default:'linear'}
  \item{nGroup}{number of group for quantile normalisation default:20}
  \item{histoneMark}{histone mark of interest. default:"ChIP-seq signal"}
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
L. Polit
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

pathToRPKMfile="/.../FPKM_values_CLBBER_CLBMA_SJNB12.txt"
sampleName=”neuroblastoma”
pathToFiles="/.../pathToFiles2.txt"
output="/.../folderTest/"

ChIPIN_normalize(RPKM=pathToRPKMfile, raw_read_count=NULL, pathToFileGenesNotMoving=NULL, sample_name=sampleName, output_dir=output, organism="hg19", fileWithPaths=pathToFiles, beforeRegionStartLength, afterRegionStartLength, regionBodyLength, binSize, expression_plot=FALSE, compute_stat=TRUE, percentage=10, typeNorm="quantile", nGroup=20)
}
