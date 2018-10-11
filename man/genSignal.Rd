\name{genSignal}
\alias{genSignal}
\title{
    Generate the Vector of SNP Effect Sizes
}
\description{
    Generate a vector of SNP effect sizes.
    Supports various signal strength distributions.
}
\usage{
genSignal(NSignalSnps, NTotalSNPs, Heritability, SignalDistr = "Same")
}
\arguments{
    \item{NSignalSnps}{
        Number of causal SNPs.
    }
    \item{SignalDistr}{
        Distribution of effect sizes. Has to be "Same", "Normal", or "Uniform".
    }
    \item{NTotalSNPs}{
        Total number of SNPs.
    }
    \item{Heritability}{
        SNP heritability, a value between 0 and 1.
    }
}
\details{
    The function generates a signal effect size vector
    that can be used in \code{\link{genSignal}}.
}
\value{
    Return a vector of length \code{NTotalSNPs}.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
%\seealso{
%    See \code{\link[base]{order}}.
%}
\examples{
genSignal(
    NSignalSnps = 10,
    NTotalSNPs = 20,
    Heritability = 0.1,
    SignalDistr = "Same")
}
