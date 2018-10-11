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
    sxy = rnorm( n = Nsnps, 
                 mean = signal * syy,
                 sd = sqrt( (1 - signal^2) * syy / df)
               );
    sxx = (1 - signal^2) / df * 
          rchisq( n = Nsnps, df = df, ncp = df * syy * signal^2 / (1 - signal^2));
                
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

