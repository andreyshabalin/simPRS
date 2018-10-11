\name{rConfInt}
\alias{rConfInt}
\title{
    Calculate Confidence Interval for Sample Correlation
}
\description{
    Calculate confidence interval for
    sample correlation using Fisher's transformation.
}
\usage{
rConfInt(r, N, alpha = 0.05)
}
\arguments{
    \item{r}{
        A vector of sample correlations. \cr
        Can be output of \code{\link{prsInf}}.
    }
    \item{N}{
        Sample size.
    }
    \item{alpha}{
        Confidence level. By default it is 5\%.
    }
}
\details{
    The function uses Fisher's transformation
    to convert sample correlations \code{r} into a z-score.
    A confidence interval is then constructed for the z-score.
    The inverse transformation is used to construct the 
    confidence interval for the sample correlations \code{r}.
}
\value{
    The function returns a data frame with 
    lower and upper bounds of the confidence intervals.
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
}
