require("survival")
require("survminer")
library("tidyverse")
library("rms")
library("patchwork")
lname <- load('./data/STAD_TCGA.RData')
lname
miseq <- stad.mirnaseq.bytype$primary.tumor %>% t() %>% as.data.frame()
dim(miseq)
rownames(miseq) <- rownames(miseq) %>% str_sub(1,12)
seq <- stad.rnaseq.bytype$primary.tumor %>% t() %>% as.data.frame()
dim(seq)
rownames(seq) <- rownames(seq) %>% str_sub(1,12)
mi_clinical <- clinical[rownames(miseq),]
seq_clinical <- clinical[rownames(seq),]
# save(seq,miseq,mi_clinical,seq_clinical,file='./mi_m_clinical.RData')

#lname <- load('../data/mi_m_clinical.RData')
#lname

# mi ----------------------------------------------------------------------

mi_clinical <- as.data.frame(mi_clinical)
mi_clinical$mir_205 <- miseq[rownames(mi_clinical),'hsa-mir-205']
mi_clinical$time <- apply(mi_clinical[,c(4,5)],1,function(x){ifelse(is.na(x[1]),x[2],x[1])}) %>% as.numeric

mi_clinical$time <- (mi_clinical$time)/30 %>% round(.,1)
mi_clinical[,c(2,3,20)] <- apply(mi_clinical[,c(2,3,20)],2,as.numeric)
mi_clinical$pathologyTstage <- mi_clinical$pathologyTstage %>% 
        str_sub(1,2) %>%
        toupper %>%
        factor(., 
               levels=c('T1','T2','T3','T4'))
mi_clinical$pathologyNstage <- mi_clinical$pathologyNstage %>% 
        as.character %>% 
        str_sub(1,2) %>% 
        toupper %>%
        factor(., 
               levels=c('N0','N1','N2','N3','NX'))
mi_clinical$pathologyMstage <- mi_clinical$pathologyMstage %>% 
        as.character %>% 
        str_sub(1,2) %>%
        toupper %>%
        factor(.,
               levels=c('M0','M1','MX'))
mi_clinical$pathologicstage <- mi_clinical$pathologicstage %>% 
        as.character %>% 
        str_replace_all(string=.,
                        pattern='stage ', 
                        '') %>% 
        str_replace_all(string=.,
                        pattern='[abc]',
                        '') %>%
        toupper %>% 
        factor(x=., 
               levels=c('I','II','III','IV'))
# m -----------------------------------------------------------------------
seq_clinical <- as.data.frame(seq_clinical)
seq_clinical$time <- apply(seq_clinical[,c(4,5)],1,function(x){ifelse(is.na(x[1]),x[2],x[1])}) 
seq_clinical[,c(2,3,19)] <- apply(seq_clinical[,c(2,3,19)],2,as.numeric)
seq_clinical$pathologyTstage <- seq_clinical$pathologyTstage %>% str_sub(1,2)

# 1 -----------------------------------------------------------------------
medians <- median(mi_clinical$mir_205)
mi_clinical$risk <- ifelse(mi_clinical$mir_205>medians, 'H','L')
fit<-survfit(Surv(time,vitalstatus)~risk,data=mi_clinical)
print(fit)

ggsurvplot(fit, palette='jco',pval=TRUE,surv.median.line='hv',cumevents=TRUE)


# therapy -----------------------------------------------------------------------
fit_therapy<-survfit(Surv(time,vitalstatus)~radiationtherapy,data=mi_clinical)
print(fit_therapy)
ggsurvplot(fit_therapy,pval =T)

# T -----------------------------------------------------------------------
fit_T<-survfit(Surv(time,vitalstatus)~pathologyTstage,data=mi_clinical)
print(fit_T)
p_sur_T <- ggsurvplot(fit_T,
                      pval =T,
                      palette=c('red','skyblue','blue','navy'),
                      surv.median.line='hv',
                      ggtheme=ggpubr::theme_pubr(),
                      pval.method=TRUE)
p_sur_T <- p_sur_T$plot+theme(legend.position='right')
# N -----------------------------------------------------------------------
fit_N<-survfit(Surv(time,vitalstatus)~pathologyNstage,data=mi_clinical)
print(fit_N)
p_sur_N <- ggsurvplot(fit_N,
                      pval =T,
                      palette='jco',
                      surv.median.line='hv',
                      ggtheme=ggpubr::theme_pubr(),
                      pval.method=TRUE)
p_sur_N <- p_sur_N$plot + theme(legend.position='right')
# M -----------------------------------------------------------------------
fit_M<-survfit(Surv(time,vitalstatus)~pathologyMstage,data=mi_clinical)
print(fit_M)
p_sur_M <- ggsurvplot(fit_M,
                      pval =T,
                      palette='jco',
                      surv.median.line='hv',
                      ggtheme=ggpubr::theme_pubr(),
                      pval.method=TRUE)
p_sur_M <- p_sur_M$plot + theme(legend.position='right')

# TNM -----------------------------------------------------------------------
fit_TNM<-survfit(Surv(time,vitalstatus)~pathologicstage,data=mi_clinical)
print(fit_TNM)
p_sur_TNM <- ggsurvplot(fit_TNM,pval =T,
                        palette='jco',
                        surv.median.line='hv',
                        ggtheme=ggpubr::theme_pubr(),
                        pval.method=TRUE)
p_sur_TNM <- p_sur_TNM$plot + theme(legend.position='right')

# Merge plots -------------------------------------------------------------
p_sur <- p_sur_TNM + {p_sur_T +
                p_sur_N +
                p_sur_M + 
                plot_layout(ncol=3)} +
        plot_layout(ncol=1)


ggsave(filename='./plots/patients_sur_TNM.pdf',
       plot=p_sur_TNM,
       device='pdf')
ggsave(filename='./plots/patients_sur_T.pdf',
       plot=p_sur_T,
       device='pdf')
ggsave(filename='./plots/patients_sur_N.pdf',
       plot=p_sur_N,
       device='pdf')
ggsave(filename='./plots/patients_sur_M.pdf',
       plot=p_sur_M,
       device='pdf')
# DE_TAR_gene ----------------------------------------------------------------------
library(pheatmap)
library(RColorBrewer)
#DE_Tar_gene <-
        # c("LDB3",
        #   "ERBB4",
        #   "SORBS1",
        #   "ZEB1",
        #   "RBPMS2",
        #   "LIMS2",
        #   "ABCD2",
        #   "RAB9B")
#DE_Tar_gene <-
        # c(
        #         "LDB3",
        #         "DSC2",
        #         "IL1RAPL1",
        #         "ERBB4"    ,
        #         "SORBS1",
        #         "CFL2",
        #         "ZEB1",
        #         "PAX9",
        #         "PROX1",
        #         "PTCHD1",
        #         "RBPMS2",
        #         "LIMS2",
        #         "NDNF",
        #         "BAMBI",
        #         "ELAVL4",
        #         "ABCD2"   ,
        #         "RAB9B",
        #         "KCNB1",
        #         "RBFOX3",
        #         "PTPRD"
        # )
        # 

gene <- readxl::read_xlsx('data/t1_VS_t234_Rna.xlsx')
DE_Tar_gene <- gene$SYMBOL
DE_Tar_gene %>% length() 

seq_clinical<- cbind(seq_clinical,seq[rownames(seq_clinical),colnames(seq) %in% DE_Tar_gene])
dim(seq_clinical)

m=0
diff_gene <- data.frame(row.names = colnames(seq_clinical)[which(colnames(seq_clinical) %in% DE_Tar_gene)])
# group <- data.frame(row.names = rownames(seq_clinical))
for(i in which(colnames(seq_clinical) %in% DE_Tar_gene)){
                n <- colnames(seq_clinical)[i]
                     print(n)
                     plot(sort(seq_clinical[, i]),
                          type = 'l',
                          main = n)
                     abline(h = quantile(seq_clinical[, i], 0.50), col = 'red')
                     text(
                             labels = '0.5',
                             col = 'red',
                             x = 200,
                             y = quantile(seq_clinical[, i], 0.80) * 1.2
                     )
                     # seq_clinical <- data.frame(seq_clinical,n=ifelse(seq_clinical[,i]>=quantile(seq_clinical[,i],0.8),'H',"L"))
                     # group <- cbind(group,n=ifelse(seq_clinical[,i]>=quantile(seq_clinical[,i],0.5),'H',"L"))
                     seq_clinical$group <-
                             ifelse(seq_clinical[, i] >= quantile(seq_clinical[, i], 0.5), 'H', "L")
                     fit_T <-
                             survfit(Surv(time, vitalstatus) ~ group, data = seq_clinical)
                     mytry <- try(fit_diff <-
                             survdiff(Surv(time, vitalstatus) ~ group, data = seq_clinical))
                     if("try-error" %in% class(mytry)){next()}
                     print(fit_diff)
                     print(i)
                     p <- 1 - pchisq(fit_diff$chisq, 1)
                     diff_gene <- rbind(diff_gene, data.frame(SYMBOL = n, pval = p))
                     if (p <= 0.05) {
                              x.inv <- try(ggsurvplot(
                                     fit_T,
                                     pval = T,
                                     title = n,
                                     surv.median.line = 'hv',
                                     risk.table = T,
                                     risk.table.y.text.col = T,
                                     ncensor.plot = T
                             ) %>% print(),
                             silent = TRUE)
                             
                             m = m + 1
                     }

                     else{
                             next()
                     }
                     
                     
}
print(m)
# [1] 959
dev.off()
save(diff_gene,file='./data/DEGs_sur.RData')
# dim(group)
# colnames(group) <- paste0(DE_Tar_gene,'_group')
# seq_clinical <- cbind(seq_clinical,group)


# sig_DE_Tar_GENE ---------------------------------------------------------
sig_DE_Tar_GENE <-
        c(      "LDB3",
                # "IL1RAPL1" 智力,
                "ERBB4",
                "SORBS1",
                "CFL2",
                "ZEB1",
                "RBPMS2",
                "LIMS2",
                "RAB9B",
                "KCNB1",
                "PTPRD"
        )
length(sig_DE_Tar_GENE)
# [1] 11

# mir-205 -----------------------------------------------------------------
summary(mi_clinical$mir_205)
plot(sort(mi_clinical$mir_205),type = 'l')
abline(h=quantile(mi_clinical$mir_205,0.50),col='red')
mi_clinical$group_by_mir_205 <- ifelse(mi_clinical$mir_205>=quantile(mi_clinical$mir_205,0.50),'H',"L")
table(mi_clinical$group_by_mir_205)
        
# H   L 
# 39 350
fit_T<-survfit(Surv(time,vitalstatus)~group_by_mir_205,data=mi_clinical)
print(fit_T)
ggsurvplot(fit_T,pval =T)

# multi GENE-------------------------------------------------------------------



# H   L 
# 59 330
fit_T1<-survfit(Surv(time,vitalstatus)~LDB3_group+SORBS1_group+ZEB1_group+RBPMS2_group+LIMS2_group+ABCD2_group,data=seq_clinical)
fit_T2<-survfit(Surv(time,vitalstatus)~LDB3_group+ZEB1_group+RBPMS2_group,data=seq_clinical)
print(fit_T1)
ggsurvplot(fit_T2,pval =T,add.all = T)



# point -------------------------------------------------------------------

load('data/DEGs_sur.RData')
diff_gene <- within(diff_gene,{
        tag <- 'black'
        tag[pval<=0.05] <- '显著'
        tag[pval>0.05] <- '不显著'
})
table(diff_gene$tag)
# 不显著   显著 
# 1498    959 
library(RColorBrewer)
p <- ggplot(diff_gene,aes(SYMBOL,pval))
p+geom_point(aes(color=tag,size=pval),alpha=0.5)+theme(axis.text.x = element_blank())+labs(x='DEG',y='P Value',color='')+scale_color_manual(values = brewer.pal(11,'RdYlBu')[c(9,11)])+geom_hline(yintercept = 0.05,color='yellow',size=1.6)+scale_size_continuous(range = c(0.1,5),limits = c(0,1))+labs(size='P Value')
ggsave(filename='plots/')
