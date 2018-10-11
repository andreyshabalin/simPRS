\name{gwasFast}
\alias{gwasFast}
\title{
    Fast Simulation of GWAS Summary Statistics
}
\description{
    Generates estimates of
    SNP effect sizes and p-values
    based on the true effect sizes and samples size.
}
\usage{
gwasFast(signal, N)
}
\arguments{
    \item{signal}{
        Vector of effect sizes for all SNPs. \cr
        Can be generated with \code{\link{genSignal}}.
    }
    \item{N}{
        GWAS sample size.
    }
}
\details{
    The function does not generate the artificial genotype and phenotype.
    Instead, finite sample distribution of the summary statistics is used
    for generation of the output.
}
\value{
    The function generates vectors of estimates of
    SNP effect sizes \code{beta} and p-values \code{pv}. \cr
    The vectors have the same length as \code{signal}.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
%\seealso{
%    See \code{\link[base]{order}}.
%}
\examples{
signal = genSignal(
            NSignalSnps = 10,
            NTotalSNPs = 20,
            heritability = 0.1,
            signalDistr = "Same")
            
gwas = gwasFast(signal = signal, N = 100)

cbind(beta = gwas$beta, pv = gwas$pv)
}
