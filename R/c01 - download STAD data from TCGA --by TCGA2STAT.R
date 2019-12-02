library(TCGA2STAT)
stad.RnaSeq2 <- getTCGA('STAD',data.type='RNASeq2',clinical=TRUE,type='RPKM')
stad.miRnaSeq2 <- getTCGA('STAD',data.type='miRNASeq',clinical=TRUE,type='rpmmm')
save(stad.RnaSeq2, file = "../data/STAD_TCGA.RData")