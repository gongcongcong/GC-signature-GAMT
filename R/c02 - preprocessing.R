library(tidyverse)
options(stringsAsFactors = FALSE)
#############################################################################
#1.input data
load('./data/STAD_TCGA.RData')
nSets = 4

DEG <- rownames(stad.rnaseq.bytype$primary.tumor)
nDEG <- length(DEG)
clinical <- as.data.frame(clinical)
setLabels =  paste0('t',1:4)
Tumor <- stad.rnaseq.bytype$primary.tumor %>% as.data.frame()
dim(Tumor)
colnames(Tumor) <- stringr::str_sub(colnames(Tumor),1,12)
myclinical <- clinical[,8] %>% as.data.frame()
myclinical$T <- stringr::str_sub(myclinical$.,1,2)
rownames(myclinical) <- rownames(clinical)
# female_index <- clinical[clinical$pathologyTstage=='female',]
# male_index <- clinical[clinical$gender=='male',]
index <- vector(mode = 'list',length =nSets )
for(i in 1:nSets) {
        index[[i]] <-
                list(index = rownames(myclinical[myclinical$T == paste0('t', i), ]))
        print(length(index[[i]]$index))
}

multiExpr = vector(mode = "list", length = nSets)

for(i in 1:nSets) {
        multiExpr[[i]] <-
                list(data = as.data.frame(apply(t(Tumor[DEG, colnames(Tumor) %in% index[[i]]$index]), 2, as.numeric)))
        
        rownames(multiExpr[[i]]$data) = colnames(Tumor[DEG, colnames(Tumor) %in% index[[i]]$index])
        
        colnames(multiExpr[[i]]$data) = rownames(Tumor[DEG, colnames(Tumor) %in% index[[i]]$index])
}
# Check that the data has the correct format 
sampleTrees = list()
for(set in 1:nSets){
        sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

######################################################################
pdf(file = "./plots/SampleClustering.pdf",
    width = 12,
    height = 12)
par(mfrow = c(2, 2))
par(mar = c(0, 4, 2, 0))
for(set in 1:nSets){
        plot(
                sampleTrees[[set]],
                main = paste("Sample clustering on all genes", setLabels[set]),
                xlab = "",
                sub = "",
                cex = 0.7
        )
}  

dev.off()


#1.c   Loading clinical trait data
traitData = STAD_RNASeq$clinical;
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = as.data.frame(traitData);
allTraits <- allTraits[rownames(allTraits) %in% colnames(Tumor),]
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)
# allTraits$race
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);

for(set in 1:nSets){
        setSamples = rownames(multiExpr[[set]]$data);
        traitRows = match(stringr::str_replace_all(setSamples,'\\.','-'), rownames(allTraits));
        Traits[[set]] = list(data = allTraits[traitRows, ]);
        rownames(Traits[[set]]$data) =setSamples;
}

ex_tumor <- stad.rnaseq.bytype$primary.tumor %>% as.data.frame
colnames(ex_tumor) <- colnames(ex_tumor) %>% str_sub(1,12)

ex_normal <- stad.rnaseq.bytype$normal %>% as.data.frame
colnames(ex_normal) <- colnames(ex_normal) %>% str_sub(1,12)
multiTraits <- Traits
# Define data set dimensions
save(multiExpr, multiTraits, clinical, ex_normal,ex_tumor,file = "./data/preprocess.RData");
