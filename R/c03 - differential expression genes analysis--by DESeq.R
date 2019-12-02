library(DESeq2)
library(tidyverse)
lname <- load('./data/STAD_TCGA.RData')
lname
lname <- load('./data/preprocess.RData')
lname
index <- c(rownames(Traits[[1]]$data),rownames(Traits[[2]]$data),rownames(Traits[[3]]$data),rownames(Traits[[4]]$data))
DE <- function(data){
        # data=stad.rnaseq.bytype
        for (i in 1:3) {
                colnames(data[[i]]) <-
                        stringr::str_sub(colnames(data[[i]]), 1, 12)
        }
        nSamples <-
                c(nrow(Traits[[1]]$data),
                  nrow(Traits[[2]]$data),nrow(Traits[[3]]$data),nrow(Traits[[4]]$data))
        group_list = factor(
                rep(paste0('t',1:4), nSamples),
                levels = paste0('t',1:4)
        )
        
        exprSet <-data[[1]][,index]
        dim(exprSet)
        colData <-
                data.frame(
                        row.names = paste(colnames(exprSet), group_list, sep = '_'),
                        group_list = group_list
                )
        dim(colData)
        colnames(exprSet) <- rownames(colData)
        rowName <- rownames(exprSet)
        exprSet <- as.data.frame(apply(exprSet, 2, as.integer))
        rownames(exprSet) <- rowName
        dim(exprSet)
        str(exprSet[1:5, 1:5])
        dds <-
                DESeqDataSetFromMatrix(
                        countData = exprSet,
                        colData = colData,
                        design = ~ group_list
                )
        
        return(dds)
        
}
lname
# dds_miRNa <- DE(data = stad.mirnaseq.bytype)
# dds2_miRNa <- DESeq(dds_miRNa,parallel=T)
dds_RNa <- DE(data=stad.rnaseq.bytype)
dds2_RNa <- DESeq(dds_RNa,parallel = T)

#load('./data/DE.RData')
t1_VS_t2_Rna <- results(dds2_RNa,contrast = c('group_list','t1','t2'),parallel = T,tidy = F)
t1_VS_t3_Rna <- results(dds2_RNa,contrast = c('group_list','t1','t3'),parallel = T,tidy = F)
t1_VS_t4_Rna <- results(dds2_RNa,contrast = c('group_list','t1','t4'),parallel = T,tidy = F)
t2_VS_t3_Rna <- results(dds2_RNa,contrast = c('group_list','t2','t3'),parallel = T,tidy = F)
t2_VS_t4_Rna <- results(dds2_RNa,contrast = c('group_list','t2','t4'),parallel = T,tidy = F)
t3_VS_t4_Rna <- results(dds2_RNa,contrast = c('group_list','t3','t4'),parallel = T,tidy = F)

# venn --------------------------------------------------------------------

A <- subset(t2_VS_t3_Rna,subset = abs(t2_VS_t3_Rna$log2FoldChange)>1&t2_VS_t3_Rna$padj<=0.05)$row %>% unlist()
B <- subset(t2_VS_t4_Rna,subset = abs(t2_VS_t4_Rna$log2FoldChange)>1&t2_VS_t4_Rna$padj<=0.05)$row %>% unlist()
C <- subset(t3_VS_t4_Rna,subset = abs(t3_VS_t4_Rna$log2FoldChange)>1&t3_VS_t4_Rna$padj<=0.05)$row %>% unlist()
length(A)
length(B)
length(C)
venn <- venn.diagram(list(t2_VS_t3=A,t2_VS_t4=B,t3_VS_t4=C),filename=NULL,col=c('red','yellow','black'),fill=c('red','yellow','black'))
grid.draw(venn)
draw.triple.venn(area1 = length(A),area2 = length(B),area3 = length(C),n12 = length(intersect(A,B)),n23 = length(intersect(B,C)),n13 =length(intersect(A,C)) ,n123 = length(intersect(A,intersect(B,C))),category = c('t2_VS_t3','t2_VS_t4','t3_VS_t4'),fill =c('red','yellow','black') )

t1_VS_t2_Rna <- as.data.frame(t1_VS_t2_Rna)
t1_VS_t3_Rna <- as.data.frame(t1_VS_t3_Rna)
t1_VS_t4_Rna <- as.data.frame(t1_VS_t4_Rna)
t2_VS_t3_Rna <- as.data.frame(t2_VS_t3_Rna)
t2_VS_t4_Rna <- as.data.frame(t2_VS_t4_Rna)
t3_VS_t4_Rna <- as.data.frame(t3_VS_t4_Rna)

t1_VS_t2_Rna$row <- rownames(t1_VS_t2_Rna)
t1_VS_t3_Rna$row <- rownames(t1_VS_t3_Rna)
t1_VS_t4_Rna$row <- rownames(t1_VS_t4_Rna)
t2_VS_t3_Rna$row <- rownames(t2_VS_t3_Rna)
t2_VS_t4_Rna$row <- rownames(t2_VS_t4_Rna)
t3_VS_t4_Rna$row <- rownames(t3_VS_t4_Rna)
save(   # dds_miRNa ,
#         dds2_miRNa,
        dds_RNa,
        dds2_RNa,
        t1_VS_t2_Rna ,
        t1_VS_t3_Rna ,
        t1_VS_t4_Rna ,
        t2_VS_t3_Rna ,
        t2_VS_t4_Rna ,
        t3_VS_t4_Rna ,
        file='../data/DE.RData'
)

# devtools::install_github("ropensci/writexl")
library(writexl)
sheets <- list("t1_VS_t2" = as.data.frame(t1_VS_t2_Rna,row.names = rownames(t1_VS_t2_Rna)),"t1_VS_t3" = as.data.frame(t1_VS_t3_Rna,row.names = rownames(t1_VS_t3_Rna)), "t1_VS_t4" = as.data.frame(t1_VS_t4_Rna,row.names = rownames(t1_VS_t4_Rna)), "t2_VS_t3" = as.data.frame(t2_VS_t3_Rna,row.names = rownames(t2_VS_t3_Rna)), "t2_VS_t4" = as.data.frame(t2_VS_t4_Rna,row.names = rownames(t2_VS_t4_Rna)), "t3_VS_t4" = as.data.frame(t3_VS_t4_Rna,row.names = rownames(t3_VS_t4_Rna))) 

write_xlsx(sheets, "./DEG/DE_RNA_20180405.xlsx")

