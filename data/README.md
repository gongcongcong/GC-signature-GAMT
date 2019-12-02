## 1.  STAD_TCGA.RData
- **description:** the data in this file, include the expression data of STAD and clinical trails of those patients
- **source:** the file was download by TCGA2STAT the from TCGA platform. 
- **code:** */R/c01 - download STAD data from TCGA --by TCGA2STAT.R*

## 2.  preprocess.RData

- **description:** to merge the clinical and the gene express pattern of each patient,  the data of STAD_TCGA.RData was preprocessed. 
- **code:** */R/c02 - preprocessing.R*

## 3.  DE_RNA_20180405.xlsx、t1_VS_t234_Rna.xlsx、DE.RData

- **description:** the differentially expression genes 
- **code:** */R/c03 - differential expression genes analysis--by DESeq.R*


## 4. DEGs_sur.RData
- **description:** Differentially expressed genes (t1 vs t2t3t4) significantly associated with survival
- **code:** */R/c04 - sur -- by survival.R*



