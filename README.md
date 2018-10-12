# simPRS: Fast Simulation of Polygenic Risk Scores

Simulation of Polygenic Risk Score calculations allows to 
(1) Select an optimal p-value threshold,
(2) Perform power analyses, and
(3) Estimate polygenicity of a phenotype.
Simulations are performed without direct simulations of the 
genotype and phenotype data. 
Instead, genome-wide association study (GWAS) summary statistics
are generated directly from their finite sample distributions.

## Installation

### Try `simPRS` online

To try `simPRS` via an online shiny interface
[click here](https://andreyshabalin.shinyapps.io/simPRS/).

### Install R package

To install `simPRS` directly from GitHub, run

```
devtools::install_github("andreyshabalin/simPRS")
```

If `devtools` package is missing, it can be installed with

```
install.packages("devtools")
```

## Sample code

### Perform a single simulation

```
# Parameters
NTotalSNPs = 10000
NSignalSnps = 100
heritability = 0.2
signalDistr = "Same"
Ntrain = 10000;
Ntest = 10000;

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
```

### Average over 1000 simulations

Let's use the parameters and `signal` defined above

```
Nsim = 1000

prsA = prsMultitest(signal = signal, N = Ntrain, Nsim = Nsim)

rci = rConfInt(r = prsA$r, N = Ntest)

prsPlot(pv = prsA$pv, r = prsA$r, rci)
```