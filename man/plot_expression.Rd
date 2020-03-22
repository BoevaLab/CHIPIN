\name{plot_expression}
\alias{plot_expression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot_expression
%%  ~~function to do ... ~~
}
\description{
Use plot_expression to profile the average ChIP-seq signal around gene TSS according to gene expression values.
}
\usage{

plot_expression(RPKM=NULL, raw_read_count=NULL, path_to_bw, organism, output_dir = ".", histone_mark = "ChIP-seq signal")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{RPKM}{path to a gene expression file (RPKM values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". If you provide the "RPKM" parameter, do not use the "raw_read_count" parameter. Default: NULL}
  \item{raw_read_count}{path to a gene expression file (raw read count values): first column should contain gene names (official gene symbol), each following column should correspond to one sample/condition. The order of values should correspond to the order of .bigWig files in "path_to_bw". Raw read count values will be transformed into RPKM values using information on exon lengths. If you provide the "raw_read_count" parameter, do not use the "RPKM" parameter. Default: NULL}
  \item{path_to_bw}{a vector containing paths to .bigWig files of the samples/conditions of interest. ! Mandatory parameter with no default value}
  \item{output_dir}{path to the output directory where one wants to store the ouput files. This directory should be created before running the function. Default: "."}
  \item{organism}{reference genome: "mm10", "mm9", "hg38" or "hg19". ! Mandatory parameter with no default value}
  \item{histone_mark}{name of the histone mark of interest; used to plot legends. Default: "ChIP-seq signal"}
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
%%  ~~who you are~~
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
outputFolder="." #change it if needed; create the corresponding output folder if it does not exists
histoneMarkName="H3K27Ac"

#run plotting function:
plot_expression(RPKM = pathToRPKMfile, path_to_bw = pathToFiles, organism="hg19", output_dir = outputFolder, histone_mark = histoneMarkName)

}

