##UniProt.ws
##availableUniprotSpecies("Drosophila melanogaster")
##  taxon ID            Species name
##1     7227 Drosophila melanogaster
##taxId(UniProt.ws) <- 7227

##Biostrings
prepareProteome <- function(UniProt.ws, fasta, species="unknown"){
    if(!missing(UniProt.ws) && class(UniProt.ws)=="UniProt.ws"){
        egs <- keys(UniProt.ws, "ENTREZ_GENE")
        cols <- c("SEQUENCE", "ID")
        proteome <- select(UniProt.ws, egs, cols, "ENTREZ_GENE")
        proteome$SEQUENCE <- gsub(" ", "", proteome$SEQUENCE,fixed=TRUE)
        proteome$LEN <- nchar(proteome$SEQUENCE)
        return(new("Proteome",
                   proteome=proteome,
                   type="UniProt",
                   species=species))
    } else {
        if(!missing(fasta)){
            if(length(fasta)==1 && class(fasta)=="character"){
                fasta <- readAAStringSet(fasta)
            }
            if(class(fasta)!="AAStringSet"){
                stop("fasta should be character or an object AAStringSet", call.=FALSE)
            }
            proteome <- data.frame(SEQUENCE=as.character(fasta),
                                   ID=names(fasta),
                                   stringsAsFactors=FALSE)
        }
        return(new("Proteome", 
                   proteome=proteome,
                   type="fasta",
                   species=species))
    }
    stop("Please check you inputs.", call.=FALSE)
}