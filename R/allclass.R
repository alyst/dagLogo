setClass("dagPeptides",
         representation(data="data.frame",
                        peptides="matrix",
                        upstreamOffset="numeric",
                        downstreamOffset="numeric",
                        type="character"),
         validity=function(object){
             re<-TRUE
             if(object@upstreamOffset < 0 || object@downstreamOffset < 0) 
                 re <- "upstreamOffset and downstreamOffset should be a integer greater than 0"
             peptides <- as.character(object@peptides)
             peptides <- peptides[(!is.na(peptides)) & (peptides!="NA")]
             if(!all(1==nchar(peptides)))
                 re <- "peptides must be a matrix with one amino acide in each position"
             re
         })

setClass("Proteome",
         representation(proteome="data.frame", type="character", species="character"),
         validity=function(object){
             re <- TRUE
             if(!object@type %in% c("fasta", "UniProt"))
                 re <- "type must be fasta or UniProt"
             if(object@type=="UniProt" && is.null(object@proteome$ENTREZ_GENE))
                 re <- "when type equals to UniProt, ENTREZ_GENE column is required for proteome"
             if(is.null(object@proteome$SEQUENCE))
                 re <- "proteome sequence is required"
             re
         })

setClass("dagBackground",
         representation(background="list", 
                        permutationSize="integer"))

setClass("testDAUresults",
         representation(group="character",
                        difference="matrix",
                        zscore="matrix",
                        pvalue="matrix",
                        background="matrix",
                        motif="matrix",
                        upstream="numeric",
                        downstream="numeric"),
         validity=function(object){
             re<-TRUE
             if(object@upstream < 0 || object@downstream < 0) 
                 re <- "upstream and downstream should be a integer greater than 0"
             
             if(ncol(object@zscore)==0 || ncol(object@difference)==0 || ncol(object@pvalue)==0)
                 re <- "slots zscore, difference and pvalue could not be empty"
             if(any(dim(object@zscore)!=dim(object@difference)) || any(dim(object@zscore)!=dim(object@pvalue)))
                 re <- "dim of slots zscore, difference and pvalue should be identical"
             re
         })
