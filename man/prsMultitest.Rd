\name{prsMultitest}
\alias{prsMultitest}
\title{
    Average Over Multiple PRS Simulations
}
\description{
    Perform a number of PRS simulations
    with \code{\link{gwasFast}} and \code{\link{prsInf}}
    and average over the results.
    Supports parallelization on multicore CPUs.
}
\usage{
prsMultitest(signal, N, Nsim, nthreads = 0, minpv = 1e-20)
}
\arguments{
    \item{signal}{
        Vector of effect sizes for all SNPs. \cr
        Can be generated with \code{\link{genSignal}}.
    }
    \item{N}{
        GWAS sample size.
    }
    \item{Nsim}{
        Number of simulations to perform.
    }
    \item{nthreads}{
        Number of parallel threads to use. \cr
        Zero (default) set it to the number of CPU cores. \cr
        Setting \code{nthreads = 1} disables multithreading.
    }
    \item{minpv}{
        The simulated PRS performance is recorded on the interval
        from \code{minpv} to 1.
    }
}
\details{
    The function uses \code{parallel} package
    for multithreaded calculations.
}
\value{
    The function returns the vector of p-value thresholds
    and corresponding average
    asymptotic correlations of
    the phenotype with PRS (\code{r}).
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
            
prsA = prsMultitest(signal = signal, N = Ntrain, Nsim = 30, nthreads = 1)

rci = rConfInt(r = prsA$r, N = Ntest)

prsPlot(pv = prsA$pv, r = prsA$r, rci)
}
