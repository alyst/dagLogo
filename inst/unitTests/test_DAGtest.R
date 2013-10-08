test_testDAU_group <- function(){
    bg <- buildBackgroundModel(seq.example, proteome=proteome.example, permutationSize=10L)
    t <- testDAU(seq.example, bg)
    checkTrue(t@zscore["K","AA0"]>100)
    checkEquals(t@pvalue["K","AA0"], 0)
    t <- testDAU(seq.example, bg, group="classic")
    checkTrue(t@zscore["positively_charged","AA0"]>50)
    checkEquals(t@pvalue["positively_charged","AA0"], 0)
}