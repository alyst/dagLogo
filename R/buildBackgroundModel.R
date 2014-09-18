buildBackgroundModel <- function(dagPeptides, 
                                 bg=c("wholeGenome", "inputSet", "nonInputSet"),
                                 model=c("any", "anchored"),
                                 targetPosition=c("any", "Nterminus", "Cterminus"),
                                 uniqueSeq=TRUE,
                                 permutationSize=30L,
                                 rand.seed=1,
                                 replacement=FALSE,
                                 proteome){
    if(missing(dagPeptides) || class(dagPeptides)!="dagPeptides"){
        stop("dagPeptides should be an object of dagPeptides.", call.=FALSE)
    }
    bg <- match.arg(bg)
    targetPosition <- match.arg(targetPosition)
    model <- match.arg(model)
    permutationSize <- as.integer(permutationSize)
    if(permutationSize<2) stop("permutationSize should be greater than 1")
    
    length <- dagPeptides@upstreamOffset + dagPeptides@downstreamOffset + 1
    ###### generate random sequences
    ## TODO, howto remove fetchSequence from background
    if(bg!="inputSet"){
        if(missing(proteome) || class(proteome)!="Proteome"){
            stop("proteome should be an object of Proteome. \n
                 Try ?prepareProteome to get help", call.=FALSE)
        }
        if(bg=="wholeGenome"){
            SequenceStr <- proteome@proteome$SEQUENCE
        }else{
            proteome.s <- proteome@proteome[!proteome@proteome$SEQUENCE %in% dagPeptides@data$peptide,]
            SequenceStr <- proteome.s$SEQUENCE
        }
    }else{
        SequenceStr <- dagPeptides@data$peptide
    }
    if(model=="anchored"){
        anchorAA <- table(dagPeptides@data$anchorAA)
        anchorAA <- anchorAA[order(anchorAA, decreasing=TRUE)]
        if(length(anchorAA)>1){
            model <- "any"
            warning("anchor amino acid is not unique. model is set to any")
            anchorAA <- paste("[",paste(names(anchorAA), collapse=""),"]", sep="")
        }else{
            anchorAA <- names(anchorAA)[1]
            pattern <- paste("([A-Z]{", dagPeptides@upstreamOffset, "}", 
                             anchorAA, "[A-Z]{", dagPeptides@downstreamOffset, "})", sep="")
        }
    }
    if(model=="any"){
        pattern <- paste("([A-Z]{", length ,"})", sep="")
    }
    if(targetPosition=="Cterminus"){
        pattern <- paste(pattern, "$", sep="")
    }else{
        if(targetPosition=="Nterminus"){
            pattern <- paste("^", pattern, sep="")
        }
    }
    matches <- gregexpr(pattern, SequenceStr)
    matches <- unlist(regmatches(SequenceStr, matches))
    set.seed(rand.seed)
    n <- nrow(dagPeptides@data)
    if(length(matches)<n) 
        stop("too less matches in background. Please try different parameters.", 
             call.=FALSE)
    background <- lapply(seq_len(permutationSize), function(p){
        s <- sample(matches, n, replace=replacement, prob=NULL)
        if(uniqueSeq){
            s <- unique(s)
        }
        do.call( rbind, strsplit( s, "", fixed = TRUE ) )
    })
    new("dagBackground", 
        background=background,
        permutationSize=permutationSize)
}