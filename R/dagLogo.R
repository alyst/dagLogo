##require pheatmap
dagHeatmap <- function(testDAUresults, type=c("diff", "zscore"), ...){
    if(missing(testDAUresults) || class(testDAUresults)!="testDAUresults"){
        stop("testDAUresults should be an object of testDAUresults\n
             Please try ?testDAU to get help.", call.=FALSE)
    }
    type <- match.arg(type)
    if(type=="diff"){
        dat <- testDAUresults@difference
    }else{
        dat <- testDAUresults@zscore
    }
    dat <- dat[order(rowSums(dat), decreasing=TRUE),]
    pheatmap(dat, ..., 
             cluster_rows=FALSE, cluster_cols=FALSE,
             scale="column")
}

colorsets <- function(colorScheme=c("null", "classic", "charge", "chemistry", "hydrophobicity")){
    colorScheme <- match.arg(colorScheme)
    auto<-c('A'='#CCFF00', 'C'='#FFFF00', 'D'='#FF0000', 'E'='#FF0066', 
            'F'='#00FF66', 'G'='#FF9900', 'H'='#0066FF', 'I'='#66FF00', 
            'K'='#6600FF', 'L'='#33FF00', 'M'='#00FF00', 'N'='#CC00FF', 
            'P'='#FFCC00', 'Q'='#FF00CC', 'R'='#0000FF', 'S'='#FF3300', 
            'T'='#FF6600', 'V'='#99FF00', 'W'='#00CCFF', 'Y'='#00FFCC')
    classic <- c("nonpolar_aliphatic"="#000000",
                    "polar_uncharged"="#00811B",
                    "aromatic"="#2000C7",
                    "positively_charged"="#800080",
                    "negatively_charged"="#FFB32C")
    chemistry <- c("hydrophobic"="#000000",
                      "polar"="#00811B",
                      "basic"="#2000C7",
                      "neutral"="#800080",
                      "acidic"='#D00001')
    hydrophobicity <- c("hydrophilic"='#000000', 
                           "neutral"='#00811B', 
                           "hydrophobic"='#2000C7')
    charge<-c("positive"="#FFB32C", "neutral"="#2000C7", "negative"="#CCCCCC")
    switch(colorScheme,
           null=auto,
           classic=classic,
           charge=charge,
           chemistry=chemistry,
           hydrophobicity=hydrophobicity,
           auto)
}

nameHash <- function(nameScheme=c("classic", "charge", "chemistry", "hydrophobicity")){
    nameScheme <- match.arg(nameScheme)
    classic <- c("nonpolar_aliphatic"="M",
                 "polar_uncharged"="U",
                 "aromatic"="A",
                 "positively_charged"="P",
                 "negatively_charged"="N")
    chemistry <- c("hydrophobic"="H",
                   "polar"="P",
                   "basic"="B",
                   "neutral"="U",
                   "acidic"='A')
    hydrophobicity <- c("hydrophilic"='I', 
                        "neutral"='U', 
                        "hydrophobic"='O')
    charge<-c("positive"="P", "neutral"="U", "negative"="N")
    switch(nameScheme,
           classic=classic,
           charge=charge,
           chemistry=chemistry,
           hydrophobicity=hydrophobicity)
}

dagLogo <- function(testDAUresults, type=c("diff", "zscore"), 
                    pvalueCutoff=0.05,
                    namehash=NULL,
                    font="Helvetica-Bold", textgp=gpar(),
                    legend=FALSE,
                    labelRelativeToAnchor=FALSE){
    if(missing(testDAUresults) || class(testDAUresults)!="testDAUresults"){
        stop("testDAUresults should be an object of testDAUresults\n
             Please try ?testDAU to get help.", call.=FALSE)
    }
    type <- match.arg(type)
    if(type=="diff"){
        dat <- testDAUresults@difference
    }else{
        dat <- testDAUresults@zscore
    }
    gscmd <- Sys.getenv("R_GSCMD")
    npos <- ncol(dat)
    ncha <- nrow(dat)
    colset <- colorsets(testDAUresults@group)
    colset <- colset[rownames(dat)]
    if(any(is.na(colset))) stop("Not every symbol has its color setting.",
                                call.=FALSE)
    if(!is.null(namehash)) {
        rownames(dat) <- namehash[rownames(dat)]
        names(colset) <- namehash[names(colset)]
    }
    rname <- rownames(dat)
    if(max(nchar(rname))>1){
        stop("Please using namehash to convert the symbols into single letters.",
             call.=FALSE)
    }
    key<-paste("x", ncha, font, paste(colset[rname], collapse=""), 
               paste(rname, collapse=""), sep="_")
    symbolsCache <- if(exists("tmp_motifStack_symbolsCache", where=".GlobalEnv")) get("tmp_motifStack_symbolsCache", pos=".GlobalEnv") else list()
    if(!is.null(symbolsCache[[key]])){
        symbols<-symbolsCache[[key]]
    } else {
        symbols<-motifStack:::coloredSymbols(ncha, font, colset[rname], rname)
        symbolsCache[[key]]<-symbols
        assign("tmp_motifStack_symbolsCache", symbolsCache, pos=".GlobalEnv")
    }
    
    dw <- 1/(npos+2)
    ##check font height to fix the number size
    plot.new()
    line1 <- as.numeric(convertUnit(unit(1, "strwidth", "W"), "npc"))
    cex <- dw/line1
    if(length(textgp)==0){
        pin <- dev.size("in")
        cex1 <- ifelse(pin[1L] > pin[2L],
                       cex * pin[2L]/pin[1L],
                       cex)
        textgp <- gpar(cex=.8 * cex1)
    }
    lwd <- cex/3
    x0 <- dw
    x1 <- 4/5 * dw
    x2 <- 6/5 * dw
    ##set ylim
    datN <- apply(dat, 2, function(.col) sum(.col[.col<0]))
    datP <- apply(dat, 2, function(.col) sum(.col[.col>0]))
    ylim <- c((as.integer(min(datN)/0.05)-1)*0.05, (as.integer(max(datP)/0.05)+1)*0.05)
    remap <- function(x){
        (ylim[1] - x)/(ylim[1] - ylim[2])/(1+dw)
    }
    reheight <- function(x){
        abs(x/(ylim[1] - ylim[2])/(1+dw))
    }
    
    ##draw axis
    x.pos <- 0
    grid.text(0, x0, remap(0)+dw/2, just=c(.5, .5), gp=textgp)
    grid.lines(x=x0, y=c(remap(0)+dw, 1), arrow=arrow(length=unit(0.02, "npc")), gp=gpar(lwd=lwd))
    grid.lines(x=x0, y=c(remap(0), 0), arrow=arrow(length=unit(0.02, "npc")), gp=gpar(lwd=lwd))
    ##draw tick
    tick <- ifelse(type=="diff", 0.1, 10^as.integer(log10(max(abs(as.numeric(dat))))))
    times <- ifelse(type=="diff", 100, 1)
    for(i in c(as.integer(min(datN)/tick):(-1), 1:as.integer(max(datP)/tick))){
        grid.lines(x=c(x2, x0), y=remap(i*tick)+dw/2, gp=gpar(lwd=lwd))
        grid.text(times*tick*i, x=x1, remap(i*tick)+dw/2, just=c(1, .5), gp=gpar(cex=lwd))
    }
    
    x.pos <- x0 + dw/2
    for(j in 1:npos){
        heights <- dat[, j]
        id <- order(heights)
        heights <- heights[testDAUresults@pvalue[, j]<pvalueCutoff]
        id <- id[id %in% which(testDAUresults@pvalue[, j]<pvalueCutoff)]
        id1 <- order(heights)
        y.pos <- remap(sum(heights[heights<0]))
        flag <- 0
        if(length(heights)>0){
            for(i in 1:length(heights)){
                h <- reheight(heights[id1[i]])
                if(heights[id1[i]]>0) flag <- flag+1
                if(flag==1) {
                    if(labelRelativeToAnchor){
                        grid.text(j-1-testDAUresults@upstream, x.pos+dw/2, y.pos+dw/2, just=c(.5, .5), gp=textgp)
                    }else{
                        grid.text(j, x.pos+dw/2, y.pos+dw/2, just=c(.5, .5), gp=textgp)
                    }
                    y.pos <- y.pos + dw
                    flag <- flag+1
                }
                if(h>0) {
                    grid.draw(grImport::pictureGrob(symbols[[id[i]]],
                                                    x.pos,y.pos,dw,h,
                                                    just=c(0,0),distort=TRUE))
                    y.pos<-y.pos+h
                }
            }
            
        }
        if(flag==0) grid.text(j, x.pos+dw/2, y.pos+dw/2, just=c(.5, .5), gp=textgp)
        x.pos<-x.pos+dw
    }
    if(legend){
        if(is.null(namehash)){
            namehash <- c("Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", "Cys"="C",
                    "Glu"="E", "Gln"="Q", "Gly"="G", "His"="H", 
                    "Ile"="I", "Leu"="L", "Lys"="K", "Met"="M", "Phe"="F", 
                    "Pro"="P", "Ser"="S", "Thr"="T", "Trp"="W", 
                    "Tyr"="Y", "Val"="V")
        }
        for(i in 1:length(namehash)){
            grid.text(namehash[i], x=0.7, y=.95-i*0.05*lwd, just=c(.5, .5), 
                      gp=gpar(col=colset[namehash[i]], cex=lwd))
            grid.text(names(namehash)[i], x=0.7+0.03*lwd, y=.95-i*0.05*lwd, 
                      just=c(0, .5), gp=gpar(col=colset[namehash[i]], cex=lwd))
        }
    }
}