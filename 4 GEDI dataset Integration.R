#===============================================================================
# *********************           Merging Datasets         *********************          
#===============================================================================


# **** Packages **** 
{
  library(GEDI)
  library(dplyr)
  library(DESeq2)
  library(tidyverse)
  library(tidyr)
  }



# **** Integration of all Dataset ****   
{
  # Character vector storing the names of directories containing individual datasets
  data_folders <- c("GSE102628", "GSE121212", "GSE140380", "GSE160501", 
                    "GSE193309", "GSE199046", "GSE230200", "GSE232127",
                    "GSE237920", "GSE256050", "GSE277961")
  
  # character vector with the data type of each dataset (RNAseq/ microarray)
  sources <- c(rep("RNAseq",11))
  
  # Reading the individual datasets using GEDI::ReadGE
  datasets <- ReadGE(data_folders, sources, path = "D:\\DRYLAB\\ATOPICDERMATITIS")
  
  # Inspecting the row identifiers of the imported datasets
  lapply(datasets, function(x){
    rownames(head(x))
  })
  
  # Attributes specifying the gene identifier type used by each dataset in the 
  # same order as the datasets listed above
  attributes <- c( rep("hgnc_symbol", 2),                   # GSE102628, GSE121212
                   rep("entrezgene_id", 5),                 # GSE140380–GSE199046
                   "ensembl_gene_id",                       # GSE232127
                   rep("entrezgene_id", 2),                 # GSE237920, GSE256050
                   "ensembl_gene_id" )                      # GSE277961
  
  
  # Integrating datasets using GEDI::GEDI using annatations from BioMart
  integrated_data <- GEDI(datasets, attributes, species = "hsapiens", BioMart = TRUE)
  

  # Constructing sample metadata for batch correction
  {
    # Collecting metadata file corresponding to each dataset
    metafiles <- list()
    for (x in data_folders) {
      folderpath <- file.path("D:/DRYLAB/ATOPICDERMATITIS", x)
      filelist <- list.files( path = folderpath,
                              pattern = "meta.*\\.rds$",
                              full.names = TRUE,
                              ignore.case = TRUE )
      metafiles[[x]] <- filelist   
      print(paste("Folder:", x, " | Meta files:", length(filelist)))
    }
    View(metafiles)
    
    # Reading the collected metadata files and retaining geo_accession no. and Skin Conditions
    metadata <- lapply(metafiles, function(file) { 
      tryCatch( 
        { 
          data <- readRDS(file) 
          data <- data[, c("geo_accession", "Condition")]  
          return(data) 
        },
        error = function(e) { 
          message(paste("Error reading file:", file, "-", e$message)) 
          NULL 
        }
      )
    })
    metadata <- bind_rows(metadata, .id = "Dataset")
    
    # Retaining metadata corresponding only to samples present in the integrated matrix
    metadata <- metadata[metadata$geo_accession %in% colnames(integrated_data), ]
    
    # Using a standard condition labels across datasets
    metadata$Condition[metadata$Condition %in% c("L", "chronic_lesion", "lesional", "atopic dermatitis", "AD", "LS")] <- "Lesional"
    metadata$Condition[metadata$Condition %in% c("NL", "non_lesional", "non-lesional", "Non-lesional", "Nonlesional")] <- "Non-lesional"
    metadata$Condition[metadata$Condition %in% c("Healthy", "HC", "healthy margin", "control", "healthy", "Normal")] <- "Healthy"
    metadata$Condition <- as.factor(metadata$Condition)
  }
  
  # Dataset specific batch assignment
  {
    # Initialising a batch vetcor
    batches <- rep(NA, ncol(integrated_data))
    
    # Filling the batch vector as per the dataset accession 
    batches <- metadata$Dataset[match(colnames(integrated_data), metadata$geo_accession)]
    
    str(batches)
    # chr [1:826] "GSE102628" "GSE102628" "GSE102628" "GSE102628" "GSE102628" ...
    # > batches
    # [1] "GSE102628" "GSE102628" "GSE102628" "GSE102628" "GSE102628" "GSE102628"
    # [7] "GSE102628" "GSE102628" "GSE102628" "GSE102628" "GSE102628" "GSE102628"
    # [13] "GSE102628" "GSE102628" "GSE102628" "GSE102628" "GSE121212" "GSE121212"
    # [19] "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212"
    # [25] "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212"
    # [31] "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212"
    # [37] "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212" "GSE121212" .....
  }
  
  # Condition specific condition assignment
  {
    # Initialising a condtion vetcor
    condition_vector <- rep(NA_character_, ncol(integrated_data))
    
    # Filling the condtion vetcor as per the dataset accession 
    condition_vector <- metadata$Condition[match(colnames(integrated_data),metadata$geo_accession)]    # Convert to factor
    condition_vector <- factor(condition_vector)
    
    str(condition_vector)
    # Factor w/ 3 levels "Healthy","Lesional",..: 2 2 2 2 2 2 2 1 1 1 ...
    #     condition_vector
    # [1] Lesional     Lesional     Lesional     Lesional     Lesional    
    # [6] Lesional     Lesional     Healthy      Healthy      Healthy     
    # [11] Healthy      Healthy      Healthy      Healthy      Healthy     
    # [16] Healthy      Lesional     Non-lesional Lesional     Non-lesional .......
    
    }
  
  # Retaining non-NA values for expression dataframe, batch vector and condition vector
  valid_idx <- !is.na(batches) & !is.na(condition_vector)
  integrated_data <- integrated_data[, valid_idx]
  batches <- batches[valid_idx]
  condition_vector <- condition_vector[valid_idx]
  
  ncol(integrated_data)
  # [1] 826
  length(batches)
  # [1] 826
  length(condition_vector)
  # [1] 826
  
  
  # Performing batch-effect correction while accounting for the biological
  # conditions and generate diagnostic visualizations using RLE boxplot
  corrected_dataset <- BatchCorrection(integrated_data, 
                                       batches, 
                                       visualize = TRUE, 
                                       status = condition_vector, 
                                       boxplot.type = "rle")
  
  MasterAD <- as.data.frame(corrected_dataset)
  saveRDS(MasterAD, "Master_Dataset.rds")
  
  
  # Verification of the corrected dataset using a tree-based approach to assess
  # whether residual structure is associated with batch or biological condition
  res <- VerifyGEDI(X = corrected_dataset, 
                    y = condition_vector, 
                    batch = batches,
                    model = "tree")
  res
  # $Feature.idx
  # [1] 10214 13252  1737  7169   856   915
  
  # $Feature.names
  # [1] "ENSG00000172243" "ENSG00000217128" "ENSG00000096060"
  # [4] "ENSG00000147601" "ENSG00000067141" "ENSG00000069275"
  
  # $TrainScores
  # GSE102628 GSE121212 GSE140380 GSE160501 GSE193309 GSE199046 GSE230200 
  # 0.6654321 0.6920981 0.6687042 0.6603774 0.6762295 0.7298050 0.6658354 
  # GSE232127 GSE237920 GSE256050 GSE277961 
  # 0.6774194 0.6687042 0.6896120 0.6927454 
  
  # $TestScores
  # GSE102628 GSE121212 GSE140380 GSE160501 GSE193309 GSE199046 GSE230200 
  # 0.7500000 0.6195652 0.6250000 0.5483871 0.6035503 0.5277778 0.7500000 
  # GSE232127 GSE237920 GSE256050 GSE277961 
  # 0.7058824 0.6250000 0.5555556 0.6097561 
  
}






