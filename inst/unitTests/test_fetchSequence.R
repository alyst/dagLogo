test_fetchSequence <- function(){
    fasta_file <- system.file("extdata", "sample_seq.fasta", package="dagLogo")
    seq_set <- readAAStringSet(fasta_file)
    names( seq_set ) <- sapply(strsplit(names(seq_set), "|", fixed=TRUE), 
                               function(chunks) chunks[[2]] ) # leave only uniprot accession
    
    dag.proteome <- prepareProteome(fasta=seq_set)
    
    dag.peptides <- fetchSequence(c('P07900', 'Q15185'), type="UniProtKB_ID", proteome=dag.proteome,
                                   anchorPos=c('Q10', 'Q84'),
                                   upstreamOffset=5, downstreamOffset=5 )
    
    checkEquals("TQTQDQPMEEE", paste(dag.peptides@peptides[1,], collapse=""))
}