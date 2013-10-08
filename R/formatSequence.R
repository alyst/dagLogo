formatSequence <- function(seq, proteome, upstreamOffset, downstreamOffset){
    if(missing(proteome) || class(proteome)!="Proteome"){
        stop("proteome should be an object of Proteome. \n
             Try ?prepareProteome to get help", call.=FALSE)
    }
    if(missing(seq)){
        stop("seq is required parameter.", call.=FALSE)
    }
    width <- unique(unlist(lapply(seq, nchar)))
    if(length(width)>1){
        stop("seq must be characters with same length", call.=FALSE)
    }
    if(missing(upstreamOffset)&&missing(downstreamOffset)){
        upstreamOffset <- floor(width/2)
        downstreamOffset <- width - upstreamOffset - 1
    }else{
        if(missing(downstreamOffset)){
            downstreamOffset <- width - upstreamOffset - 1
        }else{
            upstreamOffset <- width -downstreamOffset -1
        }
    }
    
    ## retrieve anchorAA and anchorPos
    center <- upstreamOffset + 1
    anchorAA <- unlist(lapply(seq, function(.ele) substr(.ele, center, center)))
    if(proteome@type=="UniProt") {
        type <- "entrezgene"
    } else {
        type <- ""
    }
    ## blast
    m <- do.call(rbind, lapply(seq, function(.seq) {
        .m <- regexpr(.seq, proteome@proteome$SEQUENCE)
        .id <- which(.m!=-1)
        .id <- .id[1]
        .pos <- .m[.id]
        c(.id, .pos)
    }))
    anchorPos <- m[,2]
 #   anchorAA <- anchorAA[!is.na(anchorPos)]
 #   m <- m[!is.na(anchorPos), ]
 #   anchorPos <- anchorPos[!is.na(anchorPos)]
    if(proteome@type=="UniProt") {
        IDs <- proteome@proteome[m[,1], "ENTREZ_GENE"]
    }else{
        IDs <- NA
    }
    peptide <- proteome@proteome[m[,1], "SEQUENCE"]
    anchor <- anchorAA
    dat <- data.frame(IDs=IDs, anchorAA=anchorAA, 
                      anchorPos=anchorPos, peptide=peptide, anchor=anchor)
    dat$upstream <- substr(seq, 1, upstreamOffset)
    dat$downstream <- substr(seq, upstreamOffset+2, upstreamOffset+downstreamOffset+1)
    seqchar.upstream <- do.call(rbind, lapply(dat$upstream, function(.seq){
        .seq <- c(rep("NA", upstreamOffset), 
                  unlist(lapply(1:nchar(.seq), function(i) substr(.seq, i, i))))
        .seq <- .seq[(length(.seq)-upstreamOffset+1):length(.seq)]
        .seq
    }))
    seqchar.downstream <- do.call(rbind, lapply(dat$downstream, function(.seq){
        .seq <- c(unlist(lapply(1:nchar(.seq), function(i) substr(.seq, i, i))), 
                  rep("NA", downstreamOffset))
        .seq <- .seq[1:downstreamOffset]
        .seq
    }))
    seqchar <- cbind(seqchar.upstream, as.character(dat$anchorAA), seqchar.downstream)
    new("dagPeptides", data=dat, peptides=seqchar, 
                   upstreamOffset=upstreamOffset, 
                   downstreamOffset=downstreamOffset, 
                   type=type)
    }