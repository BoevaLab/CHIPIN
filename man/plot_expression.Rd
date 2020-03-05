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

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RPKM}{dataframe with two columns: "geneName", "RPKM". If you provide RPKM parameter, put NULL for raw_read_count parameter. default:NULL}
  \item{raw_read_count}{dataframe with two columns: "geneName", "raw_read_count". Will be transformed into RPKM values using exons length information. If you provide raw_read_count parameter, put NULL for RPKM parameter. default:NULL}
  \item{fileWithPaths}{text file containing path to .bigWig files of the samples/conditions of interest. ! mandatory parameter}
  \item{output_dir}{path to output directory where one wants to store the ouput files. ! mandatory parameter}
  \item{organism}{can be mm10, mm9, hg38 or hg19. ! mandatory parameter}
  \item{histoneMark}{histone mark of interest. default:"ChIP-seq signal"
}
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
pathToFiles="/.../pathToFiles2.txt"
output="/.../folderTest/"

plot_expression(RPKM = pathToRPKMfile, raw_read_count = NULL, fileWithPaths = pathToFiles, organism=hg19, output_dir = output, histoneMark = "H3K27Ac")


}

