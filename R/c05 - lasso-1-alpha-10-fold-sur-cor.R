library(tidyverse)

# clean data --------------------------------------------------------------

lname <- load('./data/STAD_TCGA.RData')
lname <- load('./data/DEGs_sur.RData')
lname
gene <- diff_gene$SYMBOL[diff_gene$pval<=0.05]

dat <- STAD_RNASeq$merged.dat %>% as.data.frame()

dat <- dat[dat$bcr %in% rownames(clinical),c('bcr','status','OS',colnames(dat)[(colnames(dat) %in% gene)])]
dim(dat)
names(dat)

t_dat <- clinical[dat$bcr,] %>% as.data.frame()
t_dat$pathologyTstage <- str_sub(t_dat$pathologyTstage,1,2)
table(t_dat$pathologyTstage)
t_dat <- t_dat[t_dat$pathologyTstage!='tx',]
dat <- dat[dat$bcr %in% rownames(t_dat),]
dat <- dat[!is.na(dat$OS),]
rownames(dat) <- dat$bcr
# normal ------------------------------------------------------------------
mins <- apply(dat[,-c(1:3)],2,min)
maxs <- apply(dat[,-c(1:3)],2,max)
dat[,-c(1:3)] <- scale(dat[,-c(1:3)],center=mins, scale=maxs) 
# cor ---------------------------------------------------------------------
library(pheatmap)
dat_for_cor <- dat[,-c(1:3)]
phm <- cor(dat_for_cor) %>% pheatmap(color = colorRampPalette(c('green','white','red'))(100),breaks = seq(-1,1,length.out = 100),filename = './plots/corheatmap.pdf',show_rownames=F,show_colnames=F)
cluster <- cutree(phm$tree_row,2)
table(cluster) 
gene_2 <- cluster[cluster=='2'] %>% names()
write.table(gene_2,file='./plots/gene_filter_by_cor_cluster.txt')
# yangben julei  ----------------------------------------------------------
anno <- data.frame(row.names = rownames(dat),'Pathologic T stage'=t_dat[rownames(dat),8],check.names=F)
anno$`Pathologic T stage` <- anno$`Pathologic T stage` %>% toupper %>% factor(., levels=c('T1','T2','T3','T4'))
# anno$`OS > 3years or Not` <- ifelse((dat$OS >=365*3) & dat$status==0 ,'Yes', 'No')
anno$ID <- rownames(anno)
anno <- anno[order(anno$`Pathologic T stage`),]
dat <- dat[anno$ID,]
anno1 <- data.frame(row.names=anno$ID,check.names=F,'Pathologic T stage'=anno$`Pathologic T stage`)
pheatmap(dat[,gene_2],
         show_rownames =F,
         annotation_row = anno1,
         annotation_colors=list('Pathologic T stage'=c('T1'='red','T2'='skyblue','T3'='blue','T4'='navy')),
         color=colorRampPalette(colors=RColorBrewer::brewer.pal(9,'Blues'))(100),
         filename='./plots/sur_gene_cluster_heatmap2.pdf',
         show_colnames=F,
         cluster_rows=F,
         gaps_row=c(12))
# model/1 -------------------------------------------------------------------
dat_for_lasso <- dat[dat$OS>0,c('bcr','status','OS',gene_2)]
library(glmnet)
library(survival)
set.seed(20180522)
index <- sample(1:nrow(dat_for_lasso),round(0.75*nrow(dat_for_lasso)))
train <- dat_for_lasso[index,-c(1:3)]
train_y <- t_dat[dat_for_lasso[index,'bcr'],]
test <- dat_for_lasso[-index,-c(1:3)]
test_y <- t_dat[dat_for_lasso[-index,'bcr'],]
train_y$OS <- round(dat_for_lasso$OS[index]/30, 1)
test_y$OS <- round(dat_for_lasso$OS[-index]/30, 1)
fit <- cv.glmnet(x=as.matrix(train),y=Surv(time=train_y$OS, event=train_y$vitalstatus),family='cox')

plot(fit, xvar = "lambda",label=TRUE)
plot(fit, xvar = "dev", label = TRUE)

res_no_coef <- coef(fit, s = fit$lambda.min)
res_no_coef <- res_no_coef@Dimnames[[1]][which(res_no_coef !=0)]
res_no_coef
write.table(res_no_coef,'data/output/sur_cluster_10fold_1_alpha/res_sur_cluster_10fold_alhpa_1.txt')
pdf('plots/sur_cluster_1_alpha_10folds.pdf')
plot(fit)
dev.off()

# ROC ---------------------------------------------------------------------
library(survivalROC)
pre_fit <- predict(fit, newx = as.matrix(test), s = fit$lambda.min,type='response',Surv(time=test_y$OS, event=test_y$vitalstatus))
roc_pre <- survivalROC(Stime=test_y$OS, status=test_y$vitalstatus, marker = pre_fit[,1], predict.time =3*12, method="KM")
plot(roc_pre)
roc_pre_data <- data.frame(FP=roc_pre$FP,TP=roc_pre$TP)

p_ROC_sur_cluster_10fold_1_alpha <- ggplot(roc_pre_data,aes(FP,TP))+
        geom_path(color='red')+
        papaja::theme_apa()+
        geom_abline(slope=1,intercept=0,linetype='dashed')+
        ggsci::scale_color_lancet()+
        geom_text(aes(x=.1,y=.8),label=round(roc_pre$AUC,3))


ggsave('plots/p_ROC_sur_cluster_10fold_1_alpha.pdf',p_ROC_sur_cluster_10fold_1_alpha,device='pdf')
