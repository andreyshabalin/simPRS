genSignal = function(NSignalSnps, NTotalSNPs, Heritability, SignalDistr = "Same"){
    
    stopifnot( is.numeric(NSignalSnps) );
    stopifnot( length(NSignalSnps) == 1 );
    stopifnot( is.numeric(NTotalSNPs) );
    stopifnot( length(NTotalSNPs) == 1 );
    stopifnot( is.numeric(Heritability) );
    stopifnot( length(Heritability) == 1 );
    stopifnot( Heritability >= 0 );
    stopifnot( SignalDistr %in% c("Same", "Normal", "Uniform") );
    
    if( SignalDistr == 'Uniform' ){
        sig = seq(0,1, length.out = NSignalSnps+1)[-1];
    } else if( SignalDistr == 'Normal' ){
        sig = rnorm(NSignalSnps);
    } else if( SignalDistr == 'Same' ){
        sig = rep(1, NSignalSnps);
    } else {
        stop('Unknown "SignalDistr" parameter.\n',
             'Use "Same", "Normal", or "Uniform".');
    }
    sig = sqrt(abs(sig));
    sig = sig / sqrt(sum(sig^2)) * sqrt(Heritability);

    signal = c(sig, rep(0, NTotalSNPs - NSignalSnps));

    return(signal);
}