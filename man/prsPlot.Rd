\name{prsPlot}
\alias{prsPlot}
\title{
    Plot PRS performance (R-squared)
}
\description{
    Plot PRS performance (R-squared) with
    a confidence band.
}
\usage{
prsPlot(pv, r, confInt)
}
\arguments{
    \item{pv, r}{
        Vectors of p-value thresholds and 
        asymptotic correlations with phenotype,
        as produced by \code{\link{prsInf}}.
    }
    \item{confInt}{
        Upper and lower boundaries of the
        confidence band, as returned by \code{\link{rConfInt}}.
    }
}
\details{
    The plot uses reverse log10 scale for the
    p-value threshold axis (horizontal).
}
\value{
    The function returns nothing.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
%\seealso{
%    See \code{\link[base]{order}}.
%}
\examples{
NTotalSNPs = 10000
NSignalSnps = 100
heritability = 0.2
signalDistr = "Same"
Ntrain = 10000
Ntest = 3000

signal = genSignal(
            NSignalSnps = NSignalSnps,
            NTotalSNPs = NTotalSNPs,
            heritability = heritability,
            signalDistr = signalDistr)
            
gwas = gwasFast(signal = signal, N = Ntrain)

prs = prsInf(
            gwasPV = gwas$pv,
            gwasBt = gwas$beta,
            signal = signal)

rci = rConfInt(r = prs$r, N = Ntest)

prsPlot(pv = prs$pv, r = prs$r, rci)
}
