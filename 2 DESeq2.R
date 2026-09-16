#===============================================================================
#            ******** Differentially Expressed Gene Analysis ********           
#-------------------------------------------------------------------------------
#                               *** DESeq2 ***                     
#===============================================================================



# The following pipeline is the demonstration of discovery of differentially 
# expressed gene(DEG) analysis using DESeq2. The pipeline was followed for all 11 
# datasets in similarfashion. This script shows performs analysis for raw counts 
# of GSE193309. The primarygene identifier of this dataset is entrez id. The 
# supporting metadata of the datasetprovides information about the patient biopsy 
# samples namely -:skin conditions (healthy, non-lesional, lesional), 
# covariates(Anatomic region, Gender, EASI score, Scorad). The DEGs were obtained   
# both with and without covariate adjustment. 



# Packages
{
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(GEOquery)
  library(Biobase)
  library(openxlsx)
  library(tibble)
  library(DESeq2)
}



# **** Loading Data **** 
{
  GSE193309_rawcount <- read.delim("D:\\DRYLAB\\ATOPICDERMATITIS\\GSE193309\\GSE193309_raw_counts_GRCh38.p13_NCBI.tsv", header=TRUE)
  GSE193309_rawcount <- column_to_rownames(GSE193309_rawcount, "GeneID")
}



# **** Metadata ****
{
  library(GEOquery)
  clinicaldata <- getGEO(GEO="GSE193309", GSEMatrix= TRUE) # store the dataset in an object
  
  library(Biobase)
  metadata <- pData(phenoData(clinicaldata[[1]])) # extract the phenotypic data from the object using phenodata() and convert it into dataframe using pdata()
  
  meta <- metadata %>% select(c(2,12))
  names(meta)[names(meta) == "characteristics_ch1.2"] <- "Condition"
  
  meta <- meta %>% mutate(Condition = gsub("skin type:", "", meta$Condition))
  meta$Condition <- as.factor(meta$Condition)
  meta <- meta[rownames(meta) %in% colnames(GSE193309_rawcount), , drop = FALSE]
  
  # Co-variates
  {
    VarCov <- metadata[,c(2,46,47,48,49)]
    VarCov <- VarCov[VarCov$geo_accession %in% meta$geo_accession, ]
    rownames(VarCov) <- VarCov$geo_accession
    VarCov[is.na(VarCov)] <- 0
    VarCov <- VarCov %>% rename(AnatomicRegion = `anatomic_region:ch1`,
                                EASIscore = `easi:ch1`,
                                Gender = `gender:ch1`,
                                Score = `scorad:ch1`)
    
    rownames(meta) <- meta$geo_accession
    meta$geo_accession <- NULL
    VarCov$geo_accession <- NULL
    
    VarCov$EASIscore <- as.numeric(VarCov$EASIscore)
    VarCov$Score <- as.numeric(VarCov$Score)
    VarCov$AnatomicRegion <- as.factor(VarCov$AnatomicRegion)
    VarCov$Gender <- as.factor(VarCov$Gender)
    
    VarCov[is.na(VarCov)] <- 0
   }
  
  # meta + VarCov
  meta <- cbind(meta, VarCov)
  meta$geo_accession <- rownames(meta)
  meta$Score_scaled <- scale(meta$Score)
  
  saveRDS(meta, "meta.rds")
  
 }



# **** Normalisation: DESeq ****
{ 
  # Without adjustment
  {   
      # Preparing DESeqDataSet
        Dds1 <- DESeqDataSetFromMatrix(countData=GSE193309_rawcount, colData=meta, design=~Condition)
        Dds1$Condition<- relevel(Dds1$Condition, ref="HC")
        
        Dds2 <- DESeqDataSetFromMatrix(countData=GSE193309_rawcount, colData=meta, design=~Condition)
        Dds2$Condition<- relevel(Dds2$Condition, ref="NL")
        
      # Principal Component Analysis of Variance-Stabilized Expression Data and Covariate Assessment
        vsd <- vst(Dds1, blind = FALSE)
        vsd <- assay(vsd)
        plot_pca_cov(matrix_array = vsd,
                     metadata = meta,
                     covariates = c("EASIscore", "Score" ,"Gender", "AnatomicRegion","Condition"),
                     numeric_covariates = c("Score", "EASIscore"),
                     adjusted = FALSE, 
                     colours_categorical_covariate = c( "#efdf00" , "#84bd00", "#009f4d", "#00205b", "#0077c8", "#74d2e7", "#fe5000", "#da1884", "#a51890", "#ce181e", "#222"), 
                     colours_numerical_covariate = c("#004643", "#b8b8ff", "#064789"),
                     normalisation = "DESeq2")
        
      # Filter genes with raw counts >= 10 in at least 110 samples (size of the smallest group)    
        keep <- rowSums(counts(Dds1) >= 10) >= 110
        Dds1 <- Dds1[keep,]
        Dds2 <- Dds2[keep,]
        
      # Perform deseq
        Dds1 <- DESeq(Dds1) 
        Dds2 <- DESeq(Dds2)
        
      # Contrast according to the groups i.e Normal Vs Lesion and Normal Vs Non Lesional 
        resLesional <- results(Dds1, contrast = c("Condition", "LS", "HC")) %>% as.data.frame()
        resNonlesional <- results(Dds1, contrast = c("Condition", "NL", "HC")) %>% as.data.frame()
        resLesionalVsNonlesional <- results(Dds2, contrast = c("Condition", "LS", "NL")) %>% as.data.frame()
      
      # store in a list
        res <- list (  resLesional = resLesional,
                       resNonlesional = resNonlesional,
                       resLesionalVsNonlesional = resLesionalVsNonlesional)
        
      # Looping over to remove padj == NA, and add columns
        res <- lapply(res, function(x){
          x <-  x[!is.na( x$padj), ]
          x$Significant <- ifelse(x$padj < 0.05, "Significant", "Not Significant")
          x$Significant <-as.factor(x$Significant)
          x$DEG <- ifelse(x$log2FoldChange < -1, "Down" , ifelse(x$log2FoldChange > 1, "Up", "non-DEG") )
          x$DEG <-as.factor(x$DEG)
          x$Group <- paste(x$Significant, x$DEG, sep = " ")
          x$Group <-as.factor(x$Group)
          x$Gene <- rownames(x)
          x
        })
 
      # Summary of DEG analysis
        summary(res)
     
      # Saving the dataframes in rds format for meta-analysis
        saveRDS(res$resLesional, "GSE193309_Lesional.rds")
        saveRDS(res$resNonlesional, "GSE193309_Nonlesional.rds")
        saveRDS(res$resLesionalVsNonlesional, "GSE193309_LesionalVsNonlesional.rds")
      
          }
   
  # covariate adjustment
  {
      # Preparing DESeqDataSet
        DdsCov1 <- DESeqDataSetFromMatrix(countData=GSE193309_rawcount, colData=meta, design=~AnatomicRegion + Gender + Score_scaled + Condition)
        DdsCov1$Condition<- relevel(DdsCov1$Condition, ref="HC")
        
        DdsCov2 <- DESeqDataSetFromMatrix(countData=GSE193309_rawcount, colData=meta, design=~AnatomicRegion + Gender + Score_scaled + Condition)
        DdsCov2$Condition<- relevel(DdsCov2$Condition, ref="NL")
        
      # Filter genes with raw counts >= 10 in at least 110 samples (size of the smallest group)    
        keep <- rowSums(counts(Dds1_corrected) >= 10) >= 110
        DdsCov1 <- DdsCov1[keep,]
        DdsCov2 <- DdsCov2[keep,]
        
      # Extracting residuals post covariate adjustments, PCA and Covariate Assessment
        vsdCov <- vst(DdsCov1, blind = FALSE)
        deseqAdj <- limma::removeBatchEffect( assay(vsdCov),covariates = model.matrix(~ Score + AnatomicRegion + Gender, colData(vsdCov))[,-1])
        plot_pca_cov(matrix_array = deseqAdj,
                     metadata = meta,
                     covariates = c("Score" ,"Gender", "AnatomicRegion","Condition"),
                     numeric_covariates = c("Score"),
                     adjusted = TRUE, 
                     colours_categorical_covariate = c( "#efdf00" , "#84bd00", "#009f4d", "#00205b", "#0077c8","#74d2e7", "#fe5000", "#da1884", "#a51890", "#ce181e", "#222"), 
                     colours_numerical_covariate = c("#004643", "#b8b8ff", "#064789"),
                     normalisation = "DESeq2")
        
      # Perform deseq
        DdsCov1 <-DESeq(DdsCov1) 
        DdsCov2 <-DESeq(DdsCov2)
        
      # Contrast according to the groups i.e Normal Vs Lesion and Normal Vs Non Lesional 
        resCovLesional <- results(DdsCov1, contrast = c("Condition", "LS", "HC")) %>% as.data.frame()
        resCovNonlesional <- results(DdsCov1, contrast = c("Condition", "NL", "HC")) %>% as.data.frame()
        resCovLesionalVsNonlesional <- results(DdsCov2, contrast = c("Condition", "LS", "NL")) %>% as.data.frame()
    
      # store in a list
        resCov <- list (resCovLesional= resCovLesional,
                        resCovNonlesional= resCovNonlesional,
                        resCovLesionalVsNonlesional= resCovLesionalVsNonlesional )
      
      # Looping over to remove padj == NA, and add columns
        resCov <- lapply(resCov, function(x){
      x <-  x[!is.na( x$padj), ]
      x$Significant <- ifelse(x$padj < 0.05, "Significant", "Not Significant")
      x$Significant <-as.factor(x$Significant)
      x$DEG <- ifelse(x$log2FoldChange < -1, "Down" , ifelse(x$log2FoldChange > 1, "Up", "non-DEG") )
      x$DEG <-as.factor(x$DEG)
      x$Group <- paste(x$Significant, x$DEG, sep = " ")
      x$Group <-as.factor(x$Group)
      x$Gene <- rownames(x)
      x
    })
      
      # Summary of DEG analysis
        summary(resCov)
      
      # Saving the dataframes in rds format for meta-analysis  
        saveRDS(resCovLesional, "GSE193309_covLesional.rds")
        saveRDS(resCovNonlesional, "GSE193309_covNonlesional.rds")
        saveRDS(resCovLesionalVsNonlesional, "GSE193309_covLesionalVsNonlesional.rds")
      }
  
  }




