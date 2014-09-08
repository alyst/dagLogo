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
    classic <- list("nonpolar_aliphatic"=c("A", "G", "L", "M", "I", "V"),
                 "polar_uncharged"=c("C", "P", "Q", "S", "T"),
                 "aromatic"=c("F", "W", "Y"),
                 "positively_charged"=c("H", "K", "N", "R"),
                 "negatively_charged"=c("D", "E"))
    charge <- list("positive"=c("H", "K", "R"),
                "neutral"=c("A", "C", "F", "G", "I", "L", "M", "N", "P", "Q",
                            "S", "T", "V", "W", "Y"),
                "negative"=c("D", "E"))
    chemistry <- list("hydrophobic"=c("A", "F", "I", "L", "M", "P", "V", "W"),
                   "polar"=c("C", "G", "S", "T", "Y"),
                   "basic"=c("H", "K", "R"),
                   "neutral"=c("N", "Q"),
                   "acidic"=c("D", "E"))
    hydrophobicity <- list("hydrophilic"=c("D", "E", "K", "N", "Q", "R"), 
                        "neutral"=c("A", "G", "H", "P", "S", "T"), 
                        "hydrophobic"=c("C", "F", "I", "L", "M", "V", "W", "Y"))
    coln <- if(group=="null") as.character(AA) else names(get(group))
    ##convert by group
    convert <- function(x, gtype){
        d <- dim(x)
        x <- as.character(x)
        for(i in 1:length(gtype)){
            id <- x %in% gtype[[i]]
            name <- names(gtype)[i]
            x[id] <- name
        }
        matrix(x, nrow=d[1], ncol=d[2])
    }
    groupAA <- function(dat, group){
        dat <- switch(group,
                      classic=convert(dat, classic),
                      charge=convert(dat, charge),
                      hydrophobicity=convert(dat, hydrophobicity),
                      chemistry=convert(dat, chemistry),
                      null=dat
        )
        dat
    }
    bg <- lapply(bg, function(.bg) groupAA(.bg, group))
    exp <- groupAA(exp, group)
    if(ncol(exp)!=ncol(bg[[1]])){
        stop("the length of background is different from inputs", call.=FALSE)
    }
    counts <- function(mat, coln){
        num <- apply(mat, 2, function(.ele){
            cnt <- table(.ele)
            cnt <- cnt[names(cnt) %in% coln]
            total <- sum(cnt)
            percentage <- cnt/total
            percentage[coln]
        })
    }
    bg <- lapply(bg, counts, coln)
    exp <- counts(exp, coln)
    rownames(exp) <- coln
    bg <- lapply(1:ncol(exp), function(i){
        do.call(cbind, lapply(bg, function(.bg){ .bg[,i]}))
    })
    ##Z-score = (x-mu)/std
    std <- do.call(cbind, lapply(bg, function(.bg){
        apply(.bg, 1, sd, na.rm=TRUE)
    }))
    mu <- do.call(cbind, lapply(bg, function(.bg){
       apply(.bg, 1, mean, na.rm=TRUE) 
    }))
    ##difference
    exp[is.na(exp)] <- 0
    diff <- exp - mu
    diff[is.na(diff)] <- 0
    
    zscore <- diff/std
    
    rownames(diff) <- coln
    rownames(zscore) <- coln
    coln <- c()
    if(dagPeptides@upstreamOffset>0){
        coln <- paste("AA", -1*(dagPeptides@upstreamOffset:1), sep="")
    }
    coln <- c(coln, "AA0")
    if(dagPeptides@downstreamOffset>0){
        coln <- c(coln, paste("AA", 1:dagPeptides@downstreamOffset, sep=""))
    }
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
                   background=mu,
                   motif=exp,
                   upstream=dagPeptides@upstreamOffset,
                   downstream=dagPeptides@downstreamOffset)
}