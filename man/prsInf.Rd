\name{prsInf}
\alias{prsInf}
\title{
    Calculate Asymptotic Performance of 
    Poligenic Risk Score
}
\description{
    Calculate asymptotic performance of 
    poligenic risk score (PRS) at all possible p-value thresholds.
}
\usage{
prsInf(gwasBt, gwasPV, signal)
}
\arguments{
    \item{gwasBt, gwasPV}{
        Vectors of effect size estimates and p-values for all SNPs. \cr
        Can be generated with \code{\link{gwasFast}}.
    }
    \item{signal}{
        Vector of effect sizes for all SNPs. \cr
        Usually the one used to generate the GWAS summary statistics.
    }
}
\details{
    The function does not generate the artificial genotype and phenotype.
    Instead, finite sample distribution of the summary statistics is used
    for generation of the output.
}
\value{
    The function returns the vector of p-value thresholds
    (sorted vector of GWAS p-values, \code{pv})
    and corresponding asymptotic correlations of
    the phenotype with PRS (\code{r}).
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
            
gwas = gwasFast(signal = signal, N = 85)

cbind(beta = gwas$beta, pv = gwas$pv)

prs = prsInf(
            gwasPV = gwas$pv,
            gwasBt = gwas$beta,
            signal = signal)

cbind(pv = prs$pv, r = prs$r)

}
