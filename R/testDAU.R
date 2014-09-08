testDAU <- function(dagPeptides, dagBackground, 
                    group=c("null", "classic", "charge", "chemistry", "hydrophobicity")){
    if(missing(dagPeptides) || class(dagPeptides)!="dagPeptides"){
        stop("dagPeptides should be an object of dagPeptides.\n
             Please try ?fetchSequence to get help.", call.=FALSE)
    }
    if(missing(dagBackground) || class(dagBackground)!="dagBackground"){
        stop("dagBackground should be an object of dagBackground.\n
             Please try ?buildBackgroundModel to get help.", call.=FALSE)
    }
    group <- match.arg(group)
    
    exp <- dagPeptides@peptides
    bg <- dagBackground@background
    AA <- c("Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", "Cys"="C",
            "Glu"="E", "Gln"="Q", "Gly"="G", "His"="H", 
            "Ile"="I", "Leu"="L", "Lys"="K", "Met"="M", "Phe"="F", 
            "Pro"="P", "Ser"="S", "Thr"="T", "Trp"="W", 
            "Tyr"="Y", "Val"="V")
    AA_identity <- factor( AA )
    names( AA_identity ) <- AA

    aa_map <- function( class_to_aa ) {
      AA_classes <- factor( unlist( lapply( names( class_to_aa ), function( aa_class ) {
        rep.int( aa_class, length( class_to_aa[[aa_class]] ) )
      } ) ) )
      names( AA_classes ) <- unlist( class_to_aa )
      return ( AA_classes )
    }
    
    groups <- lapply( list(
    classic = list("nonpolar_aliphatic"=c("A", "G", "L", "M", "I", "V"),
                 "polar_uncharged"=c("C", "P", "Q", "S", "T"),
                 "aromatic"=c("F", "W", "Y"),
                 "positively_charged"=c("H", "K", "N", "R"),
                 "negatively_charged"=c("D", "E")),
    charge = list("positive"=c("H", "K", "R"),
                "neutral"=c("A", "C", "F", "G", "I", "L", "M", "N", "P", "Q",
                            "S", "T", "V", "W", "Y"),
                "negative"=c("D", "E")),
    chemistry = list("hydrophobic"=c("A", "F", "I", "L", "M", "P", "V", "W"),
                   "polar"=c("C", "G", "S", "T", "Y"),
                   "basic"=c("H", "K", "R"),
                   "neutral"=c("N", "Q"),
                   "acidic"=c("D", "E")),
    hydrophobicity = list("hydrophilic"=c("D", "E", "K", "N", "Q", "R"), 
                        "neutral"=c("A", "G", "H", "P", "S", "T"), 
                        "hydrophobic"=c("C", "F", "I", "L", "M", "V", "W", "Y")),
    null = AA_identity 
    ), aa_map )
    
    if(ncol(exp)!=ncol(bg[[1]])){
        stop("the length of background is different from inputs", call.=FALSE)
    }
    group_freqs <- function(x,gtype){
        df <- subset( data.frame( x = groups[[gtype]][as.character(x)],
                    col = rep( seq_len(ncol(x)), each = nrow(x) ) ), !is.na(x) )
        res <- table( df$x, df$col )
        as.matrix( as.data.frame.matrix( sweep( res, 2, colSums( res ), "/" ) ) )
    }
    bg_freqs <- lapply(bg, group_freqs, group)
    exp_freqs <- group_freqs(exp, group)
    bg_col_freqs <- lapply(seq_len(ncol(exp_freqs)), function(i){
        do.call(cbind, lapply(bg_freqs, function(.bg){ .bg[,i] }))
    })
    ##Z-score = (x-mu)/std
    bg_sd <- do.call(cbind, lapply(bg_col_freqs, function(.bg){
        apply(.bg, 1, sd, na.rm=TRUE)
    }))
    bg_mu <- do.call(cbind, lapply(bg_col_freqs, function(.bg){
       rowMeans(.bg, na.rm=TRUE)
    }))
    ##difference
    exp_freqs[is.na(exp_freqs)] <- 0
    diff <- exp_freqs - bg_mu
    diff[is.na(diff)] <- 0
    
    zscore <- diff/bg_sd

    coln <- paste("AA", seq(-dagPeptides@upstreamOffset,dagPeptides@downstreamOffset), sep="")
    if(length(coln) == ncol(diff)){
        colnames(diff) <- colnames(zscore) <- coln
    }else{
        colnames(diff) <- colnames(zscore) <- paste("AA", 1:ncol(diff), sep="")
    }
    pvalue <- 2*pnorm(-abs(zscore))  
    new("testDAUresults", group=group,
                   difference=diff,
                   zscore=zscore,
                   pvalue=pvalue,
                   background=bg_mu,
                   motif=exp,
                   upstream=dagPeptides@upstreamOffset,
                   downstream=dagPeptides@downstreamOffset)
}