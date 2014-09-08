### R code from vignette source 'dagLogo.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: fetchSequences
###################################################
library(dagLogo)
library(biomaRt)
mart <- useMart("ensembl", "dmelanogaster_gene_ensembl")
dat <- read.csv(system.file("extdata", "dagLogoTestData.csv", package="dagLogo"))
dat <- dat[1:5,] ##subset to speed sample
dat
seq <- fetchSequence(as.character(dat$entrez_geneid), 
       anchorPos=as.character(dat$NCBI_site), 
       mart=mart, 
       upstreamOffset=7, 
       downstreamOffset=7)
head(seq@peptides)


###################################################
### code chunk number 3: formatSequence
###################################################
dat <- unlist(read.delim(system.file("extdata", 
                            "grB.txt", package="dagLogo"), 
                            header=F, as.is=TRUE))
head(dat)
##prepare proteome from a fasta file
proteome <- prepareProteome(fasta=system.file("extdata", 
                                              "HUMAN.fasta",
                                              package="dagLogo"))
##prepare object of dagPeptides
seq <- formatSequence(seq=dat, proteome=proteome, 
                      upstreamOffset=14, downstreamOffset=15)


###################################################
### code chunk number 4: prepareProteome
###################################################
if(interactive()){
    library(UniProt.ws)
    taxId(UniProt.ws) <- 9606
    proteome <- prepareProteome(UniProt.ws=UniProt.ws)
}


###################################################
### code chunk number 5: prepareProteome
###################################################
bg <- buildBackgroundModel(seq, bg="wholeGenome", proteome=proteome)


###################################################
### code chunk number 6: testDAU
###################################################
t0 <- testDAU(seq, bg)
t1 <- testDAU(seq, bg, group="classic")
t2 <- testDAU(seq, bg, group="charge")
t3 <- testDAU(seq, bg, group="chemistry")
t4 <- testDAU(seq, bg, group="hydrophobicity")


###################################################
### code chunk number 7: dagHeatmap
###################################################
dagHeatmap(t0)


###################################################
### code chunk number 8: dagLogo0
###################################################
dagLogo(t0)


###################################################
### code chunk number 9: dagLogo1
###################################################
dagLogo(t1, namehash=nameHash(t1@group), legend=TRUE)


###################################################
### code chunk number 10: dagLogo2
###################################################
dagLogo(t2, namehash=nameHash(t2@group), legend=TRUE)


###################################################
### code chunk number 11: dagLogo3
###################################################
dagLogo(t3, namehash=nameHash(t3@group), legend=TRUE)


###################################################
### code chunk number 12: dagLogo4
###################################################
dagLogo(t4, namehash=nameHash(t4@group), legend=TRUE)


###################################################
### code chunk number 13: CAPmotif
###################################################
library(motifStack)
protein<-read.table(file.path(find.package("motifStack"),"extdata","cap.txt"))
protein<-t(protein[,1:20])
motif<-pcm2pfm(protein)
motif<-new("pfm", mat=motif, name="CAP", 
            color=colorset(alphabet="AA",colorScheme="chemistry"))
plot(motif)


###################################################
### code chunk number 14: CAPdagLogo
###################################################
library(Biostrings)
cap <- as.character(readAAStringSet(system.file("extdata", 
                                                "cap.fasta", 
                                                package="dagLogo")))
data(ecoli.proteome)
seq <- formatSequence(seq=cap, proteome=ecoli.proteome)
bg <- buildBackgroundModel(seq, bg="wholeGenome", 
                           proteome=ecoli.proteome, 
                           permutationSize=10L)
t0 <- testDAU(seq, bg)
dagLogo(t0)


###################################################
### code chunk number 15: CAPgroup
###################################################
t1 <- testDAU(seq, bg, group="chemistry")
dagLogo(t1, namehash=nameHash(t1@group), legend=TRUE)


###################################################
### code chunk number 16: sessionInfo
###################################################
toLatex(sessionInfo())


