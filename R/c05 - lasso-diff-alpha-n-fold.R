library(tidyverse)
library(survivalROC)
# clean data --------------------------------------------------------------

lname <- load('./data/STAD_TCGA.RData')
lname <- load('./data/DEGs_sur.RData')
lname

gene <- diff_gene$SYMBOL#[diff_gene$pval<=0.05]

dat <- STAD_RNASeq$merged.dat %>% as.data.frame()

dat <- dat[dat$bcr %in% rownames(clinical),c('bcr','status','OS',colnames(dat)[(colnames(dat) %in% gene)])]
dat <- dat[(!is.na(dat$OS))&(dat$OS!=0),]
dim(dat)
names(dat)


t_dat <- clinical[dat$bcr,] %>% as.data.frame()
t_dat$pathologyTstage <- str_sub(t_dat$pathologyTstage,1,2)
table(t_dat$pathologyTstage)
t_dat <- t_dat[t_dat$pathologyTstage!='tx',]
dat <- dat[dat$bcr %in% rownames(t_dat),]
rownames(dat) <- dat$bcr
# normal ------------------------------------------------------------------
mins <- apply(dat[,-c(1:3)],2,min)
maxs <- apply(dat[,-c(1:3)],2,max)
dat[,-c(1:3)] <- scale(dat[,-c(1:3)],center=mins, scale=maxs) 



# no operation ------------------------------------------------------------
dat_for_lasso <- dat[,as.character(gene[gene %in% colnames(dat)])]
library(glmnet)
library(survival)
set.seed(20180522)
index <- sample(1:nrow(dat_for_lasso),round(0.75*nrow(dat_for_lasso)))
train <- dat_for_lasso[index,-c(1:3)]
train_y <- t_dat[rownames(train),]
test <- dat_for_lasso[-index,-c(1:3)]
test_y <- t_dat[rownames(test),]
train_y$OS <- round(dat[rownames(train),'OS']/30, 1)
test_y$OS <- round(dat[rownames(test),'OS']/30, 1)

# 10fold, different alpha value -------------------------------------------

for (i in c(0,.2,.4,.8,1)) {
        assign(paste("cvfit", i, sep=""),
               cv.glmnet(x=as.matrix(train),y=Surv(time=train_y$OS, event=train_y$vitalstatus),family='cox',nfolds=10, alpha=i))
}
pdf('plots/no_operate_diff_alpha_10folds.pdf')
par(mfrow=c(2,3))
for (i in c(0,.2,.4,.8,1)){
        plot(get(paste('cvfit',i,sep="")),main=i)
}
dev.off()

res_no_coef0 <- coef(cvfit0, s = cvfit0$lambda.min)
res_no_coef0 <- res_no_coef0@Dimnames[[1]][which(res_no_coef0 !=0)]

res_no_coef.2 <- coef(cvfit0.2, s = cvfit0.2$lambda.min)
res_no_coef.2 <- res_no_coef.2@Dimnames[[1]][which(res_no_coef.2 !=0)]

res_no_coef.4 <- coef(cvfit0.4, s = cvfit0.4$lambda.min)
res_no_coef.4 <- res_no_coef.4@Dimnames[[1]][which(res_no_coef.4 !=0)]

res_no_coef.8 <- coef(cvfit0.8, s = cvfit0.8$lambda.min)
res_no_coef.8 <- res_no_coef.8@Dimnames[[1]][which(res_no_coef.8 !=0)]

res_no_coef1 <- coef(cvfit1, s = cvfit1$lambda.min)
res_no_coef1 <- res_no_coef1@Dimnames[[1]][which(res_no_coef1 !=0)]

print(c(res_no_coef1,res_no_coef.8,res_no_coef.4,res_no_coef.2))
write.table(res_no_coef0,'data/output/no_filter_10fold_diff_alpha/res_no_coef_alhpa_0.txt')
write.table(res_no_coef0,'data/output/no_filter_10fold_diff_alpha/res_no_coef_alhpa_0.txt')
write.table(res_no_coef.2,'data/output/no_filter_10fold_diff_alpha/res_no_coef_alhpa_.2.txt')
write.table(res_no_coef.4,'data/output/no_filter_10fold_diff_alpha/res_no_coef_alhpa_.4.txt')
write.table(res_no_coef.8,'data/output/no_filter_10fold_diff_alpha/res_no_coef_alhpa_.8.txt')
write.table(res_no_coef1,'data/output/no_filter_10fold_diff_alpha/res_no_coef_alhpa_1.txt')
# 5fold, different alpha --------------------------------------------------
for (i in c(0,.2,.4,.8,1)) {
                assign(paste("cvfit5_", i, sep=""),
                       cv.glmnet(x=as.matrix(train),y=Surv(time=train_y$OS, event=train_y$vitalstatus),family='cox',nfolds=10, alpha=i))
}
pdf('plots/no_operate_diff_alpha_5folds.pdf')
par(mfrow=c(2,3))
for (i in c(0,.2,.4,.8,1)){
        plot(get(paste('cvfit5_',i,sep="")),main=i)
}
dev.off()
res_no_coef5_0 <- coef(cvfit5_0, s = cvfit5_0$lambda.min)
res_no_coef5_0 <- res_no_coef5_0@Dimnames[[1]][which(res_no_coef5_0 !=0)]

res_no_coef5_.2 <- coef(cvfit5_0.2, s = cvfit5_0.2$lambda.min)
res_no_coef5_.2 <- res_no_coef5_.2@Dimnames[[1]][which(res_no_coef5_.2 !=0)]
res_no_coef5_.4 <- coef(cvfit5_0.4, s = cvfit5_0.4$lambda.min)
res_no_coef5_.4 <- res_no_coef5_.4@Dimnames[[1]][which(res_no_coef5_.4 !=0)]

res_no_coef5_.8 <- coef(cvfit5_0.8, s = cvfit5_0.8$lambda.min)
res_no_coef5_.8 <- res_no_coef5_.8@Dimnames[[1]][which(res_no_coef5_.8 !=0)]

res_no_coef5_1 <- coef(cvfit1, s = cvfit5_1$lambda.min)
res_no_coef5_1 <- res_no_coef5_1@Dimnames[[1]][which(res_no_coef5_1 !=0)]

write.table(res_no_coef5_0,'data/output/no_filter_5fold_diff_alpha/res_no_coef_alhpa_0.txt')
write.table(res_no_coef5_.2,'data/output/no_filter_5fold_diff_alpha/res_no_coef_alhpa_.2.txt')
write.table(res_no_coef5_.4,'data/output/no_filter_5fold_diff_alpha/res_no_coef_alhpa_.4.txt')
write.table(res_no_coef5_.8,'data/output/no_filter_5fold_diff_alpha/res_no_coef_alhpa_.8.txt')
write.table(res_no_coef5_1,'data/output/no_filter_5fold_diff_alpha/res_no_coef_alhpa_1.txt')

# fit <- glmnet(as.matrix(train),train_y$pathologyTstage,family='cox')

# AUC 10folds ---------------------------------------------------------------------
for (i in c(0,.2,.4,.8,1)) {
        assign(paste0('pre_test',i),
               predict(get(paste0('cvfit',i)), newx = as.matrix(test), s = {get(paste0('cvfit',i))}$lambda.min,type='response',Surv(time=test_y$OS, event=test_y$vitalstatus)))
}
for (i in c(0,.2,.4,.8,1)){
        assign(paste0('sur_ROC',i),
               survivalROC(Stime=test_y$OS, status=test_y$vitalstatus, marker = {get(paste0('pre_test',i))}[,1],
                           predict.time =3*12, method="KM"))
}
ROC_test_data_10nfold <- rbind(data.frame(TP =sur_ROC0$TP,FP=sur_ROC0$FP,group='sur_ROC0',label=round(sur_ROC0$AUC,3),y=.3),data.frame(TP =sur_ROC0.2$TP,FP=sur_ROC0.2$FP,group='sur_ROC0.2',label=round(sur_ROC0.2$AUC,3),y=.4))%>% 
        rbind(data.frame(TP =sur_ROC0.4$TP,FP=sur_ROC0.4$FP,group='sur_ROC0.4',label=round(sur_ROC0.4$AUC,3),y=.5))%>% 
        rbind(data.frame(TP =sur_ROC0.8$TP,FP=sur_ROC0.8$FP,group='sur_ROC0.8',label=round(sur_ROC0.8$AUC,3),y=.6))%>%
        rbind(data.frame(TP =sur_ROC1$TP,FP=sur_ROC1$FP,group='sur_ROC1',label=round(sur_ROC1$AUC,3),y=.7))

p_ROC_10fold_diff_alpha <- ggplot(ROC_test_data_10nfold,aes(FP,TP,color=group,linetype=group))+geom_path()+papaja::theme_apa()+ggsci::scale_color_lancet()+geom_text(aes(x=.8,y=y-.2,color=group,label=label))
ggsave('plots/p_ROC_10fold_diff_alpha.pdf',p_ROC_10fold_diff_alpha,device='pdf')
# AUC 5folds ---------------------------------------------------------------------
for (i in c(0,.2,.4,.8,1)) {
        assign(paste0('pre_test5_',i),
               predict(get(paste0('cvfit5_',i)), newx = as.matrix(test), s = {get(paste0('cvfit5_',i))}$lambda.min,type='response',Surv(time=test_y$OS, event=test_y$vitalstatus)))
}
for (i in c(0,.2,.4,.8,1)){
        assign(paste0('sur_ROC5_',i),
               survivalROC(Stime=test_y$OS, status=test_y$vitalstatus, marker = {get(paste0('pre_test5_',i))}[,1],
                           predict.time =3*12, method="KM"))
}
ROC_test_data_5nfold <- rbind(data.frame(TP =sur_ROC5_0$TP,FP=sur_ROC5_0$FP,group='sur_ROC5_0',label=round(sur_ROC5_0$AUC,3),y=.3),data.frame(TP =sur_ROC5_0.2$TP,FP=sur_ROC5_0.2$FP,group='sur_ROC5_0.2',label=round(sur_ROC5_0.2$AUC,3),y=.4))%>% 
        rbind(data.frame(TP =sur_ROC5_0.4$TP,FP=sur_ROC5_0.4$FP,group='sur_ROC5_0.4',label=round(sur_ROC5_0.4$AUC,3),y=.5))%>% 
        rbind(data.frame(TP =sur_ROC5_0.8$TP,FP=sur_ROC5_0.8$FP,group='sur_ROC5_0.8',label=round(sur_ROC5_0.8$AUC,3),y=.6))%>%
        rbind(data.frame(TP =sur_ROC5_1$TP,FP=sur_ROC5_1$FP,group='sur_ROC5_1',label=round(sur_ROC5_1$AUC,3),y=.7))

p_ROC_5fold_diff_alpha <- ggplot(ROC_test_data_5nfold,aes(FP,TP,color=group,linetype=group))+geom_path()+papaja::theme_apa()+ggsci::scale_color_lancet()+geom_text(aes(x=.8,y=y-.2,color=group,label=label))
ggsave('plots/p_ROC_5fold_diff_alpha.pdf',p_ROC_5fold_diff_alpha,device='pdf')
