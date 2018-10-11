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

rci = rConfInt(r = prs$r, N = 1000)

cbind(pv = prs$pv, r = prs$r, rci)

prsPlot(pv = prs$pv, r = prs$r, rci)
}
