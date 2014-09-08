##library(biomaRt)
##mart <- useMart("ensembl", "dmelanogaster_gene_ensembl")
##peptides <- fetchSequence(entrizIDs, NCBIsites, peptides, mart)
fetchSequence <- function(IDs, type="entrezgene", anchorAA=NULL, anchorPos,
                            mart, proteome, upstreamOffset, downstreamOffset){
    if(missing(mart) && missing(proteome)){
        stop("Need mart or proteome for fetching sequences.", call.=FALSE)
    }
    if(!missing(mart) && class(mart)!="Mart"){
        stop("mart should be an object of Mart", call.=FALSE)
    }
    if(!missing(proteome) && class(proteome)!="Proteome"){
        stop("proteome should be an object of Proteome. \n
                 Try ?prepareProteome to get help", call.=FALSE)
    }
    if(missing(upstreamOffset) || missing(downstreamOffset)){
        stop("Please indicate the upstreamOffset and downstreamOffset position offset anchor amino acid", call.=FALSE)
    }
    if(upstreamOffset < 0 || downstreamOffset < 0){
        stop("upstreamOffset and downstreamOffset should be a integer greater than 0", call.=FALSE)
    }
    if(downstreamOffset > 20 || upstreamOffset > 20){
        stop("upstreamOffset and downstreamOffset should be the offset of anchor amino acid, no greater than 20", call.=FALSE)
    }
    if(missing(IDs) || missing(anchorPos)){
        stop("Missing required arguments IDs or anchorPos", call.=FALSE)
    }
    if(class(IDs)!="character" || 
           !inherits(anchorPos, c("character", "numeric", "integer"))){
        stop("IDs must be characters and anchorPos should be character or number",
             call.=FALSE)
    }
    if(any(is.na(IDs)) || any(is.na(anchorPos))){
        stop("IDs or anchorPos contains NA", call.=FALSE)
    }
    if(length(IDs)!=length(anchorPos)){
        stop("length of IDs and anchorPos are not identical.", call.=FALSE)
    }
    if(class(anchorPos)=="character"){
        anchorPos <- toupper(anchorPos)
        if(any(!grepl("^[A-Z]\\d+$", anchorPos))){
            stop("anchorPos should be the amino acide followed by the position,
                 eg. K123", call.=FALSE)
        }
        anchorAA <- substr(anchorPos, 1, 1)
        anchorPos <- as.numeric(substring(anchorPos, 2))
    }
    inputs <- data.frame(ID = IDs, anchorAA = anchorAA, anchorPos = anchorPos, stringsAsFactors = FALSE)
    colnames(inputs) <- c( type, 'anchorAA', 'anchorPos' )
    ## retreive sequence
    if(!missing(mart)){
        protein <- getSequence(id=unique(as.character(IDs)), 
                               type=type,
                               seqType="peptide",
                               mart=mart)
    }else{
        id_col <- switch( type, entrezgene = 'ENTREZ_GENE', UniProtKB_ID = 'ID', 'ID' )
        if ( !id_col %in% colnames(proteome@proteome) ) {
            stop( type, " sequence identifiers are not supported by the specified proteome",
                 call.=FALSE)
        }
        protein <- proteome@proteome[proteome@proteome[,id_col] %in% unique(as.character(IDs)), c('SEQUENCE', id_col)]
        colnames(protein) <- c("peptide", type)
        protein <- protein[!is.na(protein[, type]),]
        protein$peptide <- toupper(protein$peptide)
    }
    
    dat <- merge(inputs, protein)
    dat$anchor <- unlist(mapply(substr,
                              dat$peptide, dat$anchorPos, dat$anchorPos))
    ## check that queried anchor matches the actual residue in the sequence
    if( 'anchorAA' %in% colnames(dat) ){
        anchor_mask <- dat$anchorAA == dat$anchor
        if ( any( !anchor_mask ) ) {
          warning( sum(!anchor_mask), ' sequence(s) have another AA at anchor position, excluded' )
          dat <- dat[anchor_mask,]
        }
    }
    ## extract sequences for logo
    upstreamGuard <- paste0( rep.int( "?", upstreamOffset ), collapse = '' )
    downstreamGuard <- paste0( rep.int( "?", downstreamOffset ), collapse = '' )
    peptide.guarded <- paste0( upstreamGuard, dat$peptide, downstreamGuard )
    dat$upstream <- substr(peptide.guarded, dat$anchorPos, dat$anchorPos+upstreamOffset-1)
    dat$downstream <- substr(peptide.guarded, dat$anchorPos+upstreamOffset+1, dat$anchorPos+upstreamOffset+downstreamOffset)
    # convert logo sequences into character matrix
    seqchar.upstream <- do.call(rbind, strsplit(dat$upstream, "", fixed=TRUE))
    seqchar.upstream[ seqchar.upstream == '?' ] <- NA
    seqchar.downstream <- do.call(rbind, strsplit(dat$downstream,"", fixed=TRUE))
    seqchar.downstream[ seqchar.downstream == '?' ] <- NA
    seqchar <- cbind(seqchar.upstream, dat$anchor, seqchar.downstream)
    new("dagPeptides", data=dat, peptides=seqchar, 
                   upstreamOffset=upstreamOffset, 
                   downstreamOffset=downstreamOffset, 
                   type=type)
}