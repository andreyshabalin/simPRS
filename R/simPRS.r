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