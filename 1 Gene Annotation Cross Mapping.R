#==============================================================================================================
#                                   *****  Cross-Mapping Master Table *****
#==============================================================================================================


# Packages
{
  library(readxl)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  packageVersion("org.Hs.eg.db")
  library(biomaRt)
  library(dplyr)
  library(ensembldb)
  library(AnnotationDbi)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  packageVersion("org.Hs.eg.db")
}


# Master Mapping Table
{
  ## TABLE 1------------------------------------------------------------------
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  packageVersion("org.Hs.eg.db")
  # [1] '3.20.0'

  gene_map_1 <- AnnotationDbi::select( org.Hs.eg.db,
                                      keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
                                      columns = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                                      keytype = "ENSEMBL")
  
  sum(duplicated(gene_map_1$ENSEMBL))
  # [1] 5273
  sum(duplicated(gene_map_1$ENTREZID))
  # [1] 8347
  sum(duplicated(gene_map_1$SYMBOL))
  # [1] 8347
  
  gene_map_1 <- gene_map_1 %>%
    distinct(ENSEMBL, SYMBOL, .keep_all = TRUE)
  
  # check for multipe alignment
  gene_map_1 %>% dplyr::filter(gene_map_1_clean$SYMBOL == "DDT")
  #          ,ENSEMBL ENTREZID SYMBOL
  #  1 ENSG00000099977     1652    DDT
  #  2 ENSG00000275003     1652    DDT
  
  gene_map_1%>% dplyr::filter(gene_map_1_clean$SYMBOL == "DDTL")
  #           ENSEMBL  ENTREZID SYMBOL
  # 1 ENSG00000099977 100037417   DDTL
  # 2 ENSG00000275003 100037417   DDTL
  # 3 ENSG00000099974 100037417   DDTL
  # 4 ENSG00000275758 100037417   DDTL

  ## TABLE 2------------------------------------------------------------------
  library(ensembldb)
  library(AnnotationDbi)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
  
  gene_map_2 <- AnnotationDbi::select(
    EnsDb.Hsapiens.v86,
    keys = keys(EnsDb.Hsapiens.v86, keytype = "GENEID"),
    columns = c("GENEID", "ENTREZID", "SYMBOL"),
    keytype = "GENEID")
  
  nrow(gene_map_2)
  # [1] 65059
  sum(duplicated(gene_map_2$GENEID))
  # [1] 1089
  sum(duplicated(gene_map_2$ENTREZID))
  # [1] 40006
  sum(duplicated(gene_map_2$SYMBOL))
  # [1] 8416
  
  # check for multipe alignment
  gene_map_2 %>% dplyr::filter(gene_map_2$SYMBOL == "DDT")
  #            GENEID  ENTREZID SYMBOL
  # 1 ENSG00000099977 100037417    DDT
  # 2 ENSG00000099977      1652    DDT
  # 3 ENSG00000275003 100037417    DDT
  # 4 ENSG00000275003      1652    DDT
  
  gene_map_2 %>% dplyr::filter(gene_map_2$SYMBOL == "DDTL")
  #            GENEID  ENTREZID SYMBOL
  # 1 ENSG00000099974 100037417   DDTL
  # 2 ENSG00000275758 100037417   DDTL
  
  gene_map_2 <- gene_map_2 %>%
    distinct(GENEID, ENTREZID, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(GENEID)) 
  
  nrow(gene_map_2)
  # [1] 65059
  sum(duplicated(gene_map_2$GENEID))
  # [1] 1089
  sum(duplicated(gene_map_2$ENTREZID))
  # [1] 40006
  sum(duplicated(gene_map_2$SYMBOL))
  # [1] 8416
 
  ## TABLE 3 -----------------------------------------------------------------------
  gene_map_3 <- gene_map_2 %>%
    left_join(gene_map_1%>%
                dplyr::select(ENSEMBL, SYMBOL, ENTREZID_ref = ENTREZID),
              by = c("GENEID" = "ENSEMBL", "SYMBOL")) %>%
    dplyr::mutate(ENTREZID = ifelse(!is.na(ENTREZID_ref), ENTREZID_ref, ENTREZID)) %>%
    dplyr::select(-ENTREZID_ref)
  
  # check for multipe alignment
  gene_map_3 %>% dplyr::filter(SYMBOL == "DDT")
  #            GENEID ENTREZID SYMBOL
  # 1 ENSG00000099977     1652    DDT
  # 2 ENSG00000099977     1652    DDT
  # 3 ENSG00000275003     1652    DDT
  # 4 ENSG00000275003     1652    DDT 
  
  gene_map_3 %>% dplyr::filter(SYMBOL == "DDTL")
  #            GENEID  ENTREZID SYMBOL
  # 1 ENSG00000099974 100037417   DDTL
  # 2 ENSG00000275758 100037417   DDTL
  
  # renaming the colums
  colnames(gene_map_3)[colnames(gene_map_3) == "GENEID"] <- "ENSEMBL"
  gene_map_3 <- gene_map_3[, c("ENTREZID", "ENSEMBL", "SYMBOL")]
 
  ## TABLES FROM DAVID DATABASE ----------------------------------------------------------------------------------------------------
  DAVID_1 <- read_excel("D:\\DRYLAB\\ATOPICDERMATITIS\\Mapping\\DAVIDConversion_User List_ENTREZ_GENE_ID_2026-05-05.xlsx")
  DAVID_2 <- read_excel("D:\\DRYLAB\\ATOPICDERMATITIS\\Mapping\\DAVIDConversion_User List_ENSEMBL_GENE_ID_2026-05-06.xlsx")
  DAVID_3 <- read_excel("D:\\DRYLAB\\ATOPICDERMATITIS\\Mapping\\DAVIDConversion_User List_OFFICIAL_GENE_SYMBOL_2026-05-06.xlsx")
  DAVID_4 <- read.table("D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\symbol_diff_NL mapped.txt", header = TRUE)
  
  # DAVID_1
  colnames(DAVID_1) <- c("ENTREZID", "ENSEMBL", "SYMBOL")
  DAVID_1 <- DAVID_1[, c("ENTREZID", "ENSEMBL", "SYMBOL")]
  
  # DAVID_2
  colnames(DAVID_2) <- c("ENTREZID", "ENSEMBL", "SYMBOL")
  DAVID_2 <- DAVID_2[, c("ENTREZID", "ENSEMBL", "SYMBOL")]
  
  # DAVID_3
  DAVID_3$ENSEMBL <- NA
  DAVID_3 <- DAVID_3[, c("ENTREZID", "ENSEMBL", "SYMBOL")]
  
  # DAVID_4
  colnames(DAVID_4) <- c("SYMBOL", "ENTREZID")
  DAVID_4$ENSEMBL <- NA
  DAVID_4 <- DAVID_4[, c("ENTREZID", "ENSEMBL", "SYMBOL")]
  
  
  ## TABLE 8 --------------------------------------------------------------------
  GENE_MAP <- rbind(gene_map_3, DAVID_1, DAVID_2, DAVID_3, DAVID_4) %>%
    distinct(ENSEMBL, ENTREZID, .keep_all = TRUE)
  
  # check for multipe alignment
  GENE_MAP %>% dplyr::filter(SYMBOL == "DDT")%>%
    distinct(ENSEMBL, ENTREZID, .keep_all = TRUE)
  # ENTREZID         ENSEMBL SYMBOL
  # 1     1652 ENSG00000099977    DDT
  # 2     1652 ENSG00000275003    DDT
  
  # Saving as the table for future mapping
  saveRDS(GENE_MAP, "GeneAnnotation_ensembl.v86_corrected.rds")
}

