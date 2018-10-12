genSignal = function(NSignalSnps, NTotalSNPs, heritability, signalDistr = "Same"){
    
    stopifnot( is.numeric(NSignalSnps) );
    stopifnot( length(NSignalSnps) == 1 );
    stopifnot( is.numeric(NTotalSNPs) );
    stopifnot( length(NTotalSNPs) == 1 );
    stopifnot( is.numeric(heritability) );
    stopifnot( length(heritability) == 1 );
    stopifnot( heritability >= 0 );
    stopifnot( signalDistr %in% c("Same", "Normal", "Uniform") );
    
    if( signalDistr == "Uniform" ){
        sig = seq(0,1, length.out = NSignalSnps+1)[-1];
    } else if( signalDistr == "Normal" ){
        sig = rnorm(NSignalSnps);
    } else if( signalDistr == "Same" ){
        sig = rep(1, NSignalSnps);
    } else {
        stop("Unknown \"signalDistr\" parameter.\n",
             "Use \"Same\", \"Normal\", or \"Uniform\".");
    }
    
    sig = sqrt(abs(sig));
    sig = sig / sqrt(sum(sig^2)) * sqrt(heritability);

    signal = c(sig, rep(0, NTotalSNPs - NSignalSnps));

    return(signal);
}

gwasFast = function(signal, N){
    
    stopifnot( is.numeric(signal) );
    stopifnot( is.numeric(N) );
    stopifnot( length(N) == 1 );
    
    Nsnps = length(signal);
    df = N - 1;
    syy = rchisq(n = Nsnps, df = df) / df;
    sxy = rnorm(n = Nsnps,
                mean = signal * syy,
                sd = sqrt( (1 - signal^2) * syy / df)
               );
    sxx = (1 - signal^2) / df *
          rchisq(n = Nsnps,
                 df = df, 
                 ncp = df * syy * signal^2 / (1 - signal^2));
                
    beta = sxy / sxx;
    cr = sxy / sqrt(sxx*syy);
    dfFull = N - 2;
    cor2tt = function(x){
        return( x * sqrt(dfFull / (1 - pmin(x^2,1))));
    }
    tt2pv = function(x){
        return( (pt(-abs(x), dfFull)*2) );
    }
    tt = cor2tt(cr);
    pv = tt2pv(tt);
    return(list(
            beta = beta,
            pv = pv
            # cr = cr,
            # tt = tt
        ));
}

prsInf = function(gwasBt, gwasPV, signal){
    
    stopifnot( is.numeric(gwasPV) );
    stopifnot( is.numeric(gwasBt) );
    stopifnot( is.numeric(signal) );
    
    stopifnot( length(gwasPV) == length(signal) );
    stopifnot( length(gwasBt) == length(signal) );

    
    ord = order(gwasPV);
    ord[1:100];

    ordPV = gwasPV[ord];
    ordBt = gwasBt[ord];
    ordSg = signal[ord];
    
    prsR = cumsum(ordSg * ordBt) / sqrt( cumsum( ordBt^2 ) );
    return(list(
        # ord = ord,
        pv = ordPV,
        r = prsR
        # ordLPV = rev(-log10(ordPV)),
        # ordBt = ordBt,
        # ordSg = ordSg,
        # prsR = rev(prsR)
        # prsR2 = pmax(prsR,0)^2
        ));
}

rConfInt = function(r, N, alpha = 0.05){
    
    stopifnot( is.numeric(r) );
    stopifnot( is.numeric(N) );
    stopifnot( length(N) == 1 );
    stopifnot( is.numeric(alpha) );
    stopifnot( length(alpha) == 1 );
    stopifnot( alpha <= 1 );
    stopifnot( alpha >= 0 );
    
    if( alpha > 0.5 )
        alpha = 1 - alpha;
    
    n3 = N - 3;
    z = (sqrt(n3) / 2) * log( (1 + r)/(1 - r) );
    zi = outer(z, c(-1, 1) * qnorm(alpha/2, lower.tail = FALSE), FUN = `+`);
    qi = exp(2 * zi / sqrt(n3));
    ri = (qi - 1) / (qi + 1);
    return(data.frame( 
                lower = ri[, 1], 
                upper = ri[, 2]));
}

prsPlot = function(pv, r, confInt){

    grid = -log10(pv);
    
    maxy = max(r, confInt$lower, confInt$upper)^2;
    
    # par(mgp = c(2, 1, 0));
    
    rgb = col2rgb('blue', alpha = TRUE);
    col = rgb(rgb[1]/255, rgb[2]/255, rgb[3]/255, 0.2)
    plot(
        x = grid,
        y = pmax(r,0)^2*100,
        type = 'l',
        main = "", # paste0('PRS R2, Ntest = ', param$Ntest),
        col = 'blue',
        xaxs = 'i',
        ylab = expression(paste("R"^"2"," %")),
        ylim = c(0, maxy*100), # c(0, max(ri)^2)*100,
        yaxs = 'i',
        xlab = expression(paste("-", " log"[10] * "(", italic("P"), ")")),
        xlim = c(0,max(grid)),
        las = 1);
    # mtext(expression(paste("R"^"2"," (%)")), side=2, line=2, las = 0)
    polygon(
        x = c(grid, rev(grid)), 
        y = pmax(0,c(confInt$lower, rev(confInt$upper)))^2*100,
        col = col,
        border = NA); 
    legend(
        x = 'topright', 
        legend = c(expression(paste("Asymptotic R"^"2")), 'Confidence band'),
        lwd = c(1,10),
        col = c('blue',col));
    return(invisible(NULL));
}

prsParFun = function(randseed = 0, signal, N, Nsim, minpv = 1e-20){
    
    # library(simPRS);
    
    if(randseed > 0)
        set.seed(randseed);
    
    grid = c(0.1^seq(-log10(minpv),0, by = -0.01), 1);
    
    mn = double(length(grid));
    ct = integer(length(grid));
    for( i in seq_len(Nsim) ){
        gwas = gwasFast(signal, N);
        prsI = prsInf(gwasPV = gwas$pv, gwasBt = gwas$beta, signal);
        fi = findInterval(grid, prsI$pv)
        set = which(fi > 0L);
        mn[set] = mn[set] + prsI$r[fi[set]];
        ct[set] = ct[set] + 1L;
    }
    return(list( mn = mn, ct = ct ));
}

prsMultitest = function(signal, N, Nsim, nthreads = 0, minpv = 1e-20){

    # library(parallel);
    if( nthreads == 0 )
        nthreads = detectCores(logical = TRUE);
    
    if( nthreads == 1 ){
        mnct = prsParFun(signal = signal, N = N, Nsim = Nsim, minpv = minpv)
    } else {
        cl = makeCluster(nthreads);
        # clusterExport(cl, c("gwasFast", "prsInf"));
        parlist = clusterApplyLB(cl = cl, x = 1:nthreads, fun = prsParFun,
                                 signal = signal,
                                 N = N,
                                 Nsim = ceiling(Nsim/nthreads),
                                 minpv = minpv);
        stopCluster(cl);
        mnct = list(
            mn = Reduce(`+`, lapply(parlist, `[[`, 'mn')),
            ct = Reduce(`+`, lapply(parlist, `[[`, 'ct')));
        rm(parlist, cl);
    }
    
    grid = c(0.1^seq(-log10(minpv),0, by = -0.01), 1);
    
    keep = (mnct$ct > 0.8*Nsim);
    r = mnct$mn[keep] / mnct$ct[keep];
    pv = grid[keep];
    
    return(list( pv = pv, r = r )); 
}
