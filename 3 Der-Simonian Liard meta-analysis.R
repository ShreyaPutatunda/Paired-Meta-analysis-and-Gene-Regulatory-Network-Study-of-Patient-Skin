#===============================================================================
#            ******** Differentially Expressed Gene Analysis ********           
#-------------------------------------------------------------------------------
#                    *** DerSimonian-Laird meta-analysis ***                     
#===============================================================================


# Robust DerSimonian-Laird meta-analysis for transcriptomics
# Inputs expected:
#   logFC_mat : matrix or data.frame, rows = genes, cols = studies (numeric)
#   SE_mat    : matrix or data.frame, same dim as logFC_mat, standard errors
# Returns: data.frame with per-gene meta results: pooled effect, SE, z, p, tau2,
# Q, I2, k
# dataframes with gene IDs other than entrez gene ID:-
# GSE121212 - GeneSymbol
# GSE102628 - GeneSymbol
# GSE232127 - Ensembl > GeneSymbol
# GSE277961 - Ensembl



# Functions for Der-Simonian Liard (DL) meta-analysis      
#-------------------------------------------------------------------------------
# Read all RDS objects from a specified folder
readRDS_folder <- function(folder){
  
  # Identify all RDS files in the input folder
  files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
  
  # Read each RDS file and store the resulting objects in a list
  out <- lapply(files, function(f) { readRDS(f) })
  
  # Use the corresponding filenames to label the elements of the output list
  names(out) <- basename(files)  
  
  return(out)
}


# Extracting a specified column in a harmonized dataframe from multiple data frames
extract_column_from_list <- function(df_list, column) {
  
  
  # Identify gene identifiers using the "Gene" column when available;otherwise, 
  # use the row names as gene identifiers
  get_gene_column <- function(df) {
    
    if ("Gene" %in% colnames(df)) {
      return(df$Gene)
    }
    
    if (!is.null(rownames(df))) {
      return(rownames(df))
    }
    
    stop("No gene identifiers found.")
  }
  
  
  # Process each data frame independently
  cleaned <- lapply(seq_along(df_list), function(i) {
    df <- df_list[[i]]
    
    # Construct a standardized two-column representation containing gene identifiers
    # and the requested variable
    out <- data.frame( Gene = get_gene_column(df),
                       Value = df[[column]],
                       stringsAsFactors = FALSE )
    
    # Group values by gene to identify and resolve duplicated gene identifiers
    split_list <- split(out$Value, out$Gene)
    
    # Handling multiple values for a gene in a dataframe
    collapsed <- lapply(names(split_list), function(g) {
      
      vals <- split_list[[g]]
      
      # Retaining unique non-missing values for each gene
      uniq_vals <- unique(vals[!is.na(vals)])
      
      # Case 1: all entries for the gene are missing;
      # retaining a single NA value
      if (length(uniq_vals) == 0) {
        return(data.frame(Gene = g, Value = NA))
      }
      
      # Case 2: duplicated entries contain the same non-missing value,
      # potentially accompanied by missing values;
      # collapsing these entries to a single value
      if (length(uniq_vals) == 1) {
        return(data.frame(Gene = g, Value = uniq_vals))
      }
      
      # Case 3: duplicated entries contain genuinely different non-missing values;
      # retaining all distinct values
      return(data.frame(
        Gene = rep(g, length(uniq_vals)),
        Value = uniq_vals
      ))
    })
    
    # Recombining the gene-level results into a single data frame
    out2 <- do.call(rbind, collapsed)
    
    # Renaming the extracted value column using the corresponding dataset name
    names(out2)[2] <- names(df_list)[i]
    
    return(out2)
  })
  
  
  # Merge the processed datasets by gene identifier while retaining genes 
  # present in any dataset
  merged <- Reduce(function(x, y){
    merge(x, y, by = "Gene", all = TRUE)},cleaned)
  
  return(merged)
}


# Helper: single-gene DL meta-analysis
meta_dl_single_gene <- function(beta, 
                           se, 
                           method = c("random", "fixed"), 
                           min_se = 1e-8) {
  
  # Selecting the requested meta-analysis model: "random" or "fixed"
  method <- match.arg(method)
  
  # Retaining only observations with finite effect estimates and standard errors
  # This prevents NA, NaN, Inf, or -Inf values from affecting the calculations
  ok <- is.finite(beta) & is.finite(se)
  beta <- beta[ok]
  se   <- se[ok]
  
  # Number of studies/effect estimates contributing to the meta-analysis
  k <- length(beta)
  
  # If no valid studies/effect estimates remain, return NA for all
  if (k == 0) {
    return(list(
      beta = NA_real_, 
      se = NA_real_, 
      z = NA_real_, 
      p = NA_real_,
      Q = NA_real_, 
      tau2 = NA_real_, 
      I2 = NA_real_, 
      k = 0L
    ))
  }
  
  # Replacing zero or negative standard errors with a small positive value
  # to avoid division by zero when calculating inverse-variance weights
  se[se <= 0] <- min_se
  
  # Calculating the fixed-effect inverse-variance weights
  # Studies with smaller standard errors receive greater weight
  w <- 1 / (se^2)
  print("fixed effect weights")
  print(w)
  
  # Calculating the fixed-effect pooled estimate using inverse-variance weighting
  beta_fe <- sum(w * beta) / sum(w)
  print("fixed effect pooled effect")
  print(beta_fe)
  
  # Cochran's Q statistic measures the weighted deviation of individual study 
  # effects from the fixed-effect pooled estimate
  # Calculating Cochran's Q statistic to assess between-study heterogeneity
  Q <- sum(w * (beta - beta_fe)^2)
  print("Q")
  print(Q)
  
  # Degrees of freedom for Cochran's Q
  df <- k - 1
  print("Degree of Freedom")
  print(df)
  
  # Calculating the C term used in the DerSimonian-Laird estimator of 
  # between-study variance (tau^2)
  C <- sum(w) - sum(w^2) / sum(w)
  print("C")
  print(C)
  
  # Estimating between-study variance (tau^2) using the DerSimonian-Laird(DL) method
  # Initialising tau2
  tau2 <- 0
  
  # Caluculating Raw DL estimate of between-study variance
  # Cases when tau2 is a finite value or remains zero
  if (k > 1) {
    tau2_raw <- (Q - df) / C
    print("tau2")
    print(tau2_raw)
    # Constrain tau^2 to zero when the raw estimate is negative,zero, or non-finite, 
    # since a variance cannot be negative
    tau2 <- ifelse(is.finite(tau2_raw) & tau2_raw > 0, tau2_raw, 0)
      } else {
    # Between-study variance cannot be estimated from a single study i.e k=1
    tau2 <- 0
    print("tau2")
    print(tau2)
  }
  
  # Calculating Weights as per method chosen or tau2 value
  # Using fixed-effect weights when a fixed-effect model is requested or when the 
  # estimated between-study variance is zero. When tau^2 = 0, the random-effects 
  # model reduces to the fixed-effect model.
  if (method == "fixed" || tau2 == 0) {
    # Fixed-effect inverse-variance weights
    w_star <- w
    print("w*")
    print(w_star)
  } else {
    # Random-effects weights incorporating both within-study variance (SE^2) and 
    # between-study variance (tau^2)
    w_star <- 1 / (se^2 + tau2)
    print("w*")
    print(w_star)
  } 
  
  # Calculating the final pooled effect estimate using the selected fixed-effect
  # or random-effects weights
  beta_hat <- sum(w_star * beta) / sum(w_star)
  print("beta_hat")
  print(beta_hat)
  
  # Calculating the final standard error of the pooled effect estimate
  se_hat <- sqrt(1 / sum(w_star))
  print("se_hat")
  print(se_hat)
  
  # Calculating the Z-statistic for testing whether the pooled effect differs 
  # from zero
  z <- beta_hat / se_hat
  print("z")
  print(z)
  
  # Calculating the two-sided P-value corresponding to the Z-statistic
  pval <- 2 * pnorm(-abs(z))
  print("pval")
  print(pval)
  
  # Estimating I^2, the percentage of total observed variation attributable to
  # between-study heterogeneity rather than sampling 
  # Initialisng
  I2 <- 0
  
  # Caluculating I^2
  # Cases when I^2 is a finite %  or remains zero
  if (k > 1 && Q > df) {
    # DerSimonian-Laird estimate of I^2, expressed as a percentage
    I2 <- max(0, (Q - df)/Q ) * 100
    print("I2")
    print(I2) } else {
    # Set I^2 to zero when there is only one study or when
    # Q does not exceed its degrees of freedom
    I2 <- 0
    print("I2")
    print(I2)
  }
  
  # Return the pooled effect, its uncertainty, heterogeneity statistics,
  # and the number of contributing studies/effect estimates
  return(list(
    beta = as.numeric(beta_hat),
    se   = as.numeric(se_hat),
    z    = as.numeric(z),
    p    = as.numeric(pval),
    Q    = as.numeric(Q),
    tau2 = as.numeric(tau2),
    I2   = as.numeric(I2),
    k    = as.integer(k)
  ))
}


# Vectorized wrapper: gene-matrix DL meta-analysis
meta_dl_gene_matrix <- function(logFC_mat, 
                                SE_mat, 
                                genes = NULL, 
                                method = c("random", "fixed"),
                                min_se = 1e-8, 
                                p_adjust = c("none", "BH", "bonferroni")) {
  # Selecting the requested meta-analysis model: "random" or "fixed"
  method <- match.arg(method)
  
  # Selecting the requested p adjustment method: "BH" for Benjamini-Hochberg FDR correction,
  # "bonferroni" or "none"
  p_adjust <- match.arg(p_adjust)
  
  # Converting the input effect-size and standard-error dataframes to matrices
  logFC_mat <- as.matrix(logFC_mat)
  SE_mat    <- as.matrix(SE_mat)
  
  # Check : effect size and standard-error should have same dimensions
  if (!all(dim(logFC_mat) == dim(SE_mat))) stop("logFC_mat and SE_mat must have same dimensions.")
  
  # Determining the number of genes (rows) and studies/datasets (columns)
  n_genes <- nrow(logFC_mat)
  n_stud  <- ncol(logFC_mat)
  
  # Assigning the gene names
  # If gene names are not provided, use row names from the logFC matrix.
  # If row names are also unavailable, generate default gene names.
  if (is.null(genes)) {
    genes <- if (!is.null(rownames(logFC_mat))) {
      rownames(logFC_mat)
    } else {
      paste0("Gene", seq_len(n_genes))
    }
  }
  
  # Ensure that the number of gene names matches the number of genes
  if (length(genes) != n_genes) {
    stop("genes length must match number of rows in logFC_mat.")
  }
  
  # Allocating result vectors
  res_beta <- numeric(n_genes)   # pooled effect estimate
  res_se   <- numeric(n_genes)   # standard error of pooled effect
  res_z    <- numeric(n_genes)   # Z-statistic
  res_p    <- numeric(n_genes)   # unadjusted P-value
  res_Q    <- numeric(n_genes)   # Cochran's Q statistic
  res_tau2 <- numeric(n_genes)   # between-study variance
  res_I2   <- numeric(n_genes)   # percentage heterogeneity
  res_k    <- integer(n_genes)   # number of contributing studies
  
  
  # Performing the meta-analysis independently on each gene.(Looping DL over genes)
  for (i in seq_len(n_genes)) {
    
    # Extracting the effect estimates (logFC) and standard errors for the 
    # current gene across all studies/datasets
    beta_i <- logFC_mat[i, ]
    se_i   <- SE_mat[i, ]
    
    # Applying the single-gene DerSimonian-Laird meta-analysis functionusing the
    # selected fixed- or random-effects model
    out <- meta_dl_single( beta = beta_i, 
                           se = se_i, 
                           method = method, 
                           min_se = min_se)
    
    # Storing the meta-analysis results for the current gene in the corresponding 
    # pre-allocated result vectors
    res_beta[i] <- out$beta
    res_se[i]   <- out$se
    res_z[i]    <- out$z
    res_p[i]    <- out$p
    res_Q[i]    <- out$Q
    res_tau2[i] <- out$tau2
    res_I2[i]   <- out$I2
    res_k[i]    <- out$k
  }
  
  
  # Adjust the gene-level P-values for multiple testing. No adjustment is  
  # applied when p_adjust = "none".
  if (p_adjust == "none") {
    adjp <- res_p
    } else {
    # Applying the selected multiple-testing correction across all genes
    adjp <- p.adjust(  res_p, 
                       method = ifelse(p_adjust == "BH", "BH", "bonferroni"))
  }
  
  
  # Combining the meta-analysis results into a single data frame,
  res_df <- data.frame( gene = genes,
                        k = res_k,
                        beta = res_beta,
                        se = res_se,
                        z = res_z,
                        p = res_p,
                        p_adj = adjp,
                        Q = res_Q,
                        tau2 = res_tau2,
                        I2 = res_I2,
                        stringsAsFactors = FALSE )
  
  # Removing row names 
  rownames(res_df) <- NULL
  
  # Return the complete gene-level meta-analysis results
  return(res_df)
}




# 1. Cross Mapping to Entrez gene id
#-------------------------------------------------------------------------------
GeneAnnotation_ensemble.v86_corrected <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\GeneAnnotation_ensembl.v86_corrected.rds")

# Lesional vs Healthy
{
  # Loading DESeq2 analysis result for every dataset from respective folder
  folder_L <- "D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\Lesional"
  Lesional <- readRDS_folder(folder_L)
  
  # Coverage of Mapping table
  length(setdiff(rownames(Lesional_new$GSE277961_Lesional_corrected.rds), GeneAnnotation_ensemble.v86_corrected$ENSEMBL))
  # [1] 331
  length(setdiff(rownames(Lesional_new$GSE232127_Lesional.rds), GeneAnnotation_ensemble.v86_corrected$ENSEMBL))
  # [1] 109
  length(setdiff(rownames(Lesional_new$GSE121212_Lesional.rds), GeneAnnotation_ensemble.v86_corrected$SYMBOL))
  # [1] 1137
  length(setdiff(rownames(Lesional_new$GSE102628_Lesional.rds), GeneAnnotation_ensemble.v86_corrected$SYMBOL))
  # [1] 853  
  length(setdiff(rownames(Lesional_new$GS140380_Lesional_corrected.rds), GeneAnnotation_ensemble.v86_corrected$ENTREZID))
  # [1] 635
  
  # Mapping to entrez id
  # GSE102628
  Lesional$GSE102628_Lesional.rds <- Lesional$GSE102628_Lesional.rds  %>% 
    mutate( gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(Gene, GeneAnnotation_ensemble.v86_corrected$SYMBOL)],
            gene_id = ifelse(is.na(gene_id), Gene, gene_id)) %>%
    rename(SYMBOL = Gene)
  
  # GSE121212
  Lesional$GSE121212_Lesional.rds <- Lesional$GSE121212_Lesional.rds  %>%
    mutate(SYMBOL = rownames(.)) %>% 
    mutate( gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(SYMBOL, GeneAnnotation_ensemble.v86_corrected$SYMBOL)],
            gene_id = ifelse(is.na(gene_id), SYMBOL, gene_id))
  
  # GSE232127
  Lesional$GSE232127_Lesional.rds <-  Lesional$GSE232127_Lesional.rds %>% 
                                      rownames_to_column("GeneID") %>%
                                      separate(GeneID, into = c("ENSEMBL", "SYMBOL"), sep = ">")
  rownames(Lesional$GSE232127_Lesional.rds) <- Lesional$GSE232127_Lesional.rds$ENSEMBL
  Lesional$GSE232127_Lesional.rds <- Lesional$GSE232127_Lesional.rds %>%
    mutate(gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(ENSEMBL, GeneAnnotation_ensemble.v86_corrected$ENSEMBL)],
           gene_id = ifelse(is.na(gene_id), ENSEMBL, gene_id))
 
  # GSE277961
  Lesional$GSE277961_Lesional_corrected.rds <- Lesional$GSE277961_Lesional_corrected.rds %>% 
    mutate(ENSEMBL = rownames(.)) %>% 
    mutate( gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(ENSEMBL, GeneAnnotation_ensemble.v86_corrected$ENSEMBL)],
            gene_id = ifelse(is.na(gene_id), ENSEMBL, gene_id))
  
  # Standardising the gene identifier column before meta-analysis
  Lesional <- lapply(Lesional, function(df) {
    
    # Adding Gene column only if absent
    if (!"Gene" %in% colnames(df)) {
      
      # Using gene_id if present
      if ("gene_id" %in% colnames(df)) {
        df$Gene <- df$gene_id
        
        # Otherwise using row names
      } else {
        df$Gene <- rownames(df)
      }
    }
    
    return(df)
  })
}

# Non-lesional vs Healthy
{ 
  # Loading DESeq2 analysis result for every dataset from respective folder
  folder_NL <- "D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\NonLesional\\New"
  NonLesional <- readRDS_folder(folder_NL)
  
  # Coverage of Mapping table
  length(setdiff(rownames(NonLesional$GSE277961_Nonlesional_corrected.rds), GeneAnnotation_ensemble.v86_corrected$ENSEMBL))
  # [1] 331
  length(setdiff(rownames(NonLesional$GSE232127_Nonlesional.rds), GeneAnnotation_ensemble.v86_corrected$ENSEMBL))
  # [1] 109
  length(setdiff(rownames(NonLesional$GSE121212_Nonlesional.rds), GeneAnnotation_ensemble.v86_corrected$SYMBOL))
  # [1] 778
  length(setdiff(rownames(NonLesional$GSE121212_Nonlesional.rds), GeneAnnotation_ensemble.v86_corrected$SYMBOL))
  # [1] 772
  length(setdiff(rownames(NonLesional$GSE230200_Nonlesional_corrected.rds), GeneAnnotation_ensemble.v86_corrected$ENTREZID))
  # [1] 679
  
  # Mapping to entrez id
  # GSE121212
  NonLesional$GSE121212_Nonlesional.rds <- NonLesional$GSE121212_Nonlesional.rds  %>%
    mutate(SYMBOL = rownames(.)) %>% 
    mutate( gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(SYMBOL, GeneAnnotation_ensemble.v86_corrected$SYMBOL)],
            gene_id = ifelse(is.na(gene_id), SYMBOL, gene_id))
  
  # GSE277961
  NonLesional$GSE277961_Nonlesional_corrected.rds <-  NonLesional$GSE277961_Nonlesional_corrected.rds %>% 
    mutate(ENSEMBL = rownames(.)) %>% 
    mutate( gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(ENSEMBL, GeneAnnotation_ensemble.v86_corrected$ENSEMBL)],
            gene_id = ifelse(is.na(gene_id), ENSEMBL, gene_id))
  
  # GSE232127
  NonLesional$GSE232127_Nonlesional.rds <- NonLesional$GSE232127_Nonlesional.rds%>% 
    rownames_to_column("GeneID") %>%
    separate(GeneID, into = c("ENSEMBL", "SYMBOL"), sep = ">")
  rownames(NonLesional$GSE232127_Nonlesional.rds) <- NonLesional$GSE232127_Nonlesional.rds$ENSEMBL
  
  NonLesional$GSE232127_Nonlesional.rds <- NonLesional$GSE232127_Nonlesional.rds %>%
    mutate(gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(ENSEMBL, GeneAnnotation_ensemble.v86_corrected$ENSEMBL)],
           gene_id = ifelse(is.na(gene_id), ENSEMBL, gene_id))
  
  # Standardising the gene identifier column before meta-analysis
  NonLesional <- lapply(NonLesional, function(df) {
    
    # Adding Gene column only if absent
    if (!"Gene" %in% colnames(df)) {
      
      # Using gene_id if present
      if ("gene_id" %in% colnames(df)) {
        df$Gene <- df$gene_id
        
        # Otherwise using row names
      } else {
        df$Gene <- rownames(df)
      }
    }
    
    return(df)
  })
  
}

# Lesional vs Non-lesional
{
  # Loading DESeq2 analysis result for every dataset from respective folder
  folder_LNL <- "D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\LesionalVsNon-Lesional\\New"
  LesionalVsNonLesional <- readRDS_folder(folder_new_LNL)
  
  # Coverage of Mapping table
  length(setdiff(rownames(LesionalVsNonLesional$GSE277961_LesionalVsNonlesional_corrected.rds), GeneAnnotation_ensemble.v86_corrected$ENSEMBL))
  # [1] 331
  length(setdiff(rownames(LesionalVsNonLesional$GSE232127_LesionalVsNonlesional.rds), GeneAnnotation_ensemble.v86_corrected$ENSEMBL))
  # [1] 109
  length(setdiff(rownames(LesionalVsNonLesional$GSE121212_LesionalVsNonlesional.rds), GeneAnnotation_ensemble.v86_corrected$SYMBOL))
  # [1] 772
  length(setdiff(rownames(LesionalVsNonLesional$GSE230200_LesionalVsNonlesional_corrected.rds), GeneAnnotation_ensemble.v86_corrected$ENTREZID))
  # [1] 679
  
  # Mapping to entrez id
  # GSE121212
  LesionalVsNonLesional$GSE121212_LesionalVsNonlesional.rds <- LesionalVsNonLesional$GSE121212_LesionalVsNonlesional.rds  %>%
    mutate(SYMBOL = rownames(.)) %>% 
    mutate( gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(SYMBOL, GeneAnnotation_ensemble.v86_corrected$SYMBOL)],
            gene_id = ifelse(is.na(gene_id), SYMBOL, gene_id))
  
  # GSE232127
  LesionalVsNonLesional$GSE232127_LesionalVsNonlesional.rds <- LesionalVsNonLesional$GSE232127_LesionalVsNonlesional.rds %>% 
    rownames_to_column("GeneID") %>%
    separate(GeneID, into = c("ENSEMBL", "SYMBOL"), sep = ">")
  rownames(LesionalVsNonLesional$GSE232127_LesionalVsNonlesional.rds) <- LesionalVsNonLesional$GSE232127_LesionalVsNonlesional.rds$ENSEMBL
  LesionalVsNonLesional$GSE232127_LesionalVsNonlesional.rds <- LesionalVsNonLesional$GSE232127_LesionalVsNonlesional.rds %>%
    mutate(gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(ENSEMBL, GeneAnnotation_ensemble.v86_corrected$ENSEMBL)],
           gene_id = ifelse(is.na(gene_id), ENSEMBL, gene_id))
  
  # GSE277961
  LesionalVsNonLesional$GSE277961_LesionalVsNonlesional_corrected.rds <-  LesionalVsNonLesional$GSE277961_LesionalVsNonlesional_corrected.rds %>% 
    mutate(ENSEMBL = rownames(.)) %>% 
    mutate( gene_id = GeneAnnotation_ensemble.v86_corrected$ENTREZID[match(ENSEMBL, GeneAnnotation_ensemble.v86_corrected$ENSEMBL)],
            gene_id = ifelse(is.na(gene_id), ENSEMBL, gene_id))
  
  # Standardising the gene identifier column before meta-analysis
  LesionalVsNonLesional <- lapply(LesionalVsNonLesional, function(df) {
    
    # Adding Gene column only if absent
    if (!"Gene" %in% colnames(df)) {
      
      # Using gene_id if present
      if ("gene_id" %in% colnames(df)) {
        df$Gene <- df$gene_id
        
        # Otherwise using row names
      } else {
        df$Gene <- rownames(df)
      }
    }
    
    return(df)
  })  
}



# 2. Gene-wise Meta-Analysis across 11 datasets
#-------------------------------------------------------------------------------
# Lesional Vs Healthy 
{
  # Extracting log2 fold changes and their corresponding standard errors from 
  # each Lesional vs Healthy DESeq2 result
  Lesional_LFC <- extract_column_from_list(Lesional, "log2FoldChange")
  Lesional_SE <- extract_column_from_list(Lesional, "lfcSE")

  # Define the gene identifiers used to label the meta-analysis results
  genes <- Lesional_LFC$Gene
  
  # Convert the extracted log2 fold changes and standard errors into matrices
  # with genes as rows and studies/datasets as columns
  logFC_mat <- as.matrix(Lesional_LFC[, -1])
  storage.mode(logFC_mat) <- "numeric"
  SE_mat    <- as.matrix(Lesional_SE[, -1])
  storage.mode(SE_mat) <- "numeric"
  
  # Perform gene-wise random-effects DerSimonian-Laird meta-analysis across datasets
  # P-values are adjusted for multiple testing using the Benjamini-Hochberg procedure
  meta_L_deg <- meta_dl_matrix(logFC_mat = logFC_mat,
                               SE_mat = SE_mat,
                               genes = genes,
                               method = "random",
                               min_se = 1e-8,
                               p_adjust = "BH")
  meta_L_deg <- meta_L_deg %>% filter( k >= 10,               # Significant DEGs
                                       p_adj < 0.05,
                                       beta >= 1 | beta <= -1,
                                       !is.na(z))
  
  # Top differentially expressed genes with no heterogeneity
  Low_het_topL <- meta_L_deg %>% filter(I2 ==0) %>% arrange(desc(beta))
  
  # Annotating with gene symbols
  meta_L_deg <- meta_L_deg %>% mutate(GeneSymbol = GeneAnnotation_ensemble.v86_corrected$SYMBOL[match(gene, GeneAnnotation_ensemble.v86_corrected$ENTREZID)])
  rownames(meta_L_deg) <- meta_L_deg$GeneSymbol
  
  # Inspecting duplicates
  dup_L <- meta_L_deg$gene[duplicated(meta_L_deg$gene)==TRUE]
  meta_L_deg %>% filter(gene %in% dup_L)
  #   gene  k     beta        se        z            p        p_adj        Q      tau2       I2 GeneSymbol
  # 1 9997 12 1.498835 0.1694652 8.844501 9.193832e-19 6.995936e-17 47.77253 0.1965809 76.97422       SCO2
  # 2 9997 12 1.516823 0.1598257 9.490480 2.299733e-21 2.394865e-19 47.83036 0.1831416 77.00205       SCO2
  
  write.xlsx( meta_L_deg, "meta-analysis_L_DEG.xlsx", overwrite = TRUE)
  saveRDS(meta_L_deg,"meta-analysis_L_DEG.rds" )
  write.xlsx( Low_het_topL_new, "Low heterogeniety top L.xlsx", overwrite = TRUE)
}


# Non-lesional Vs Healthy 
{
  # Extracting log2 fold changes and their corresponding standard errors from 
  # each Non-esional vs Healthy DESeq2 result
  NonLesional_LFC <- extract_column_from_list(NonLesional, "log2FoldChange")
  NonLesional_SE <- extract_column_from_list(NonLesional, "lfcSE")
  
  # Define the gene identifiers used to label the meta-analysis results
  genes <- NonLesional_LFC$Gene
  
  # Convert the extracted log2 fold changes and standard errors into matrices
  # with genes as rows and studies/datasets as columns
  logFC_mat <- as.matrix(NonLesional_LFC[, -1])
  storage.mode(logFC_mat) <- "numeric"
  SE_mat    <- as.matrix(NonLesional_SE[, -1])
  storage.mode(SE_mat) <- "numeric"
  
  # Perform gene-wise random-effects DerSimonian-Laird meta-analysis across datasets
  # P-values are adjusted for multiple testing using the Benjamini-Hochberg procedure
  meta_all_NL_deg <- meta_dl_matrix(logFC_mat = logFC_mat,
                                         SE_mat = SE_mat,
                                         genes = genes,
                                         method = "random",
                                         min_se = 1e-8,
                                         p_adjust = "BH")
   meta_NL_deg <- meta_all_NL_deg %>% filter( k >= 4,         # Significant DEGs
                                             p_adj < 0.05,
                                             beta >= 1 | beta <= -1,
                                             !is.na(z))
  
  # Top differentially expressed genes with no heterogeneity
  Low_het_topNL_new <- meta_NL_deg %>% filter(I2 ==0) %>% arrange(desc(beta))
  
  # Annotating with gene symbols
  meta_NL_deg <- meta_NL_deg %>% mutate(GeneSymbol = GeneAnnotation_ensemble.v86_corrected$SYMBOL[match(gene, GeneAnnotation_ensemble.v86_corrected$ENTREZID)])
  rownames(meta_NL_deg) <- meta_NL_deg$GeneSymbol
  
  write.xlsx( meta_NL_deg, "meta-analysis_NL_DEG.xlsx", overwrite = TRUE)
  saveRDS(meta_NL_deg,"meta-analysis_NL_DEG.rds" )
  write.xlsx( Low_het_topNL, "Low heterogeniety top NL.xlsx", overwrite = TRUE)
}


# Lesional Vs Non-lesional 
{
  # Extracting log2 fold changes and their corresponding standard errors from 
  # each Lesional vs Healthy DESeq2 result
  LesionalVsNonLesional_LFC <- extract_column_from_list(LesionalVsNonLesional_new, "log2FoldChange")
  LesionalVsNonLesional_SE <- extract_column_from_list(LesionalVsNonLesional_new, "lfcSE")
  
  # Convert the extracted log2 fold changes and standard errors into matrices
  # with genes as rows and studies/datasets as columns
  logFC_mat <- as.matrix(LesionalVsNonLesional_LFC[, -1])
  storage.mode(logFC_mat) <- "numeric"
  SE_mat    <- as.matrix(LesionalVsNonLesional_SE[, -1])
  storage.mode(SE_mat) <- "numeric"
  
  
  # Define the gene identifiers used to label the meta-analysis results
  genes <- LesionalVsNonLesional_LFC$Gene
  

  # Perform gene-wise random-effects DerSimonian-Laird meta-analysis across datasets
  # P-values are adjusted for multiple testing using the Benjamini-Hochberg procedure
  meta_LNL_deg <- meta_dl_matrix(logFC_mat = logFC_mat,
                                     SE_mat = SE_mat,
                                     genes = genes,
                                     method = "random",
                                     min_se = 1e-8,
                                     p_adjust = "BH")
  meta_LNL_deg <- meta_LNL_deg %>% filter( k >= 4,            # Significant DEGs
                                           p_adj < 0.05,
                                           beta >= 1 | beta <= -1,
                                           !is.na(z))
  
  # Top differentially expressed genes with no heterogeneity
  Low_het_topLNL_new <- meta_LNL_deg_new %>% filter(I2 ==0) %>% arrange(desc(beta))
  
  # Annotating with gene symbols
  meta_LNL_deg <- meta_LNL_deg %>%
    mutate(GeneSymbol = GeneAnnotation_ensemble.v86_corrected$SYMBOL[match(gene, GeneAnnotation_ensemble.v86_corrected$ENTREZID)]) %>%
    mutate(GeneSymbol = ifelse(!is.na(GeneSymbol), GeneSymbol, gene))
  rownames(meta_LNL_deg_new) <- meta_LNL_deg_new$gene
  
  # Inspecting duplicates
  dup_LNL <- meta_LNL_deg[duplicated(meta_LNL_deg$GeneSymbol), ]
  meta_LNL_deg %>% filter(GeneSymbol == "SCO2")
  # gene k     beta        se        z          p      p_adj         Q tau2 I2 GeneSymbol
  # 9997 7 1.107784 0.4442076 2.493844 0.01263681 0.02268357 0.7807428    0  0       SCO2
  # 9997 7 1.180868 0.4792667 2.463905 0.01374326 0.02446798 0.9291731    0  0       SCO2
  meta_LNL_deg %>% filter(GeneSymbol == "FAM95B1")
  #      gene k      beta        se         z            p       p_adj            Q      tau2  I2 GeneSymbol
  # 100133036 4 -1.031133 0.3118189 -3.306833 9.435719e-04 0.001927832 1.166772e+16 0.3889242 100    FAM95B1
  # 105379252 5 -1.029172 0.2543012 -4.047059 5.186509e-05 0.000116329 1.293382e+16 0.3233454 100    FAM95B1
  meta_LNL_deg %>% filter(GeneSymbol == "RP11-475O6.1")
  #        gene k      beta        se         z            p        p_adj            Q      tau2  I2   GeneSymbol
  # 1 101927560 4 -1.053453 0.2800675 -3.761426 1.689476e-04 3.664209e-04 9.412534e+15 0.3137511 100 RP11-475O6.1
  # 2 101927587 6 -1.044905 0.2341518 -4.462508 8.100578e-06 1.898778e-05 1.315701e+16 0.3289252 100 RP11-475O6.1
  
  
  write.xlsx( meta_LNL_deg_new, "meta-analysis_LNL_DEG_new.xlsx", overwrite = TRUE)
  saveRDS(meta_LNL_deg_new,"meta-analysis_LNL_DEG_new.rds" )
  write.xlsx( Low_het_topLNL_new, "Low heterogeniety top LNL new.xlsx", overwrite = TRUE)
}



