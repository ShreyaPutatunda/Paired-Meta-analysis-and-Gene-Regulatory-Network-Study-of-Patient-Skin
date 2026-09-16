# ==============================================================================
# ********************************** FUNCTIONS *********************************
# ==============================================================================



# 1.plot_pca_cov
#-------------------------------------------------------------------------------
# Performs principal component analysis (PCA) on an expression matrix and 
# visualizes sample distribution according to selected categorical and/or  
# continuous metadata covariates. The function can be used to compare sample  
# structure before or after covariate adjustment.
# Returns a ggplot2 PCA object displaying PC1 versus PC2, faceted by the  
# selected covariates.
# matrix_array = dataframe; metadata = dataframe; covariates = vector pass as string
# numeric_covariates = vector; adjustment_status? = before/ after;
# colours_categorical_covariate = vector in the order of values;
# colours_numerical_covariate = vector in order of low, mid, high 

plot_pca_cov <- function(matrix_array,
                         metadata , 
                         covariates, 
                         numeric_covariates = NULL, 
                         adjusted = FALSE, 
                         colours_categorical_covariate = NULL, 
                         colours_numerical_covariate = NULL,
                         normalisation){
  if(inherits(matrix_array, c("matrix", "array", "data.frame"))){
    
    any(grepl("^GSM|^Sample", colnames(matrix_array)))
    if(any(grepl("^GSM|^Sample", colnames(matrix_array)))){
      
      library(dplyr)
      library(tidyr)
      library(tidyverse)
      library(ggplot2)
      library(ggpubr)
      library(ggnewscale)
      library(uwot)
      
      
      # Analysis
      PCA <- prcomp(t(matrix_array))
      PCAscores <- PCA$x %>% as.data.frame() 
      
      # Modifying the score dataframe
      { 
        rownames(PCAscores) <- rownames(t(matrix_array)) 
        PCAscores$samples <- rownames(PCAscores)            
        PCAscores <- full_join(PCAscores, metadata, by = c("samples" = "geo_accession"))
        
        # reshaping score matrix to long format for selected covariates
        meta_long <- PCAscores %>%
          mutate(across(covariates, as.character)) %>%
          pivot_longer( cols = covariates,
                        names_to = "Covariate",
                        values_to = "Value") %>%
          mutate( CovType = ifelse(Covariate %in% c(numeric_covariates), "continuous", "categorical"))
        all_of(covariates)
      }
      
      # Handle missing numeric covariates
      if (is.null(numeric_covariates)) { numeric_covariates <- character(0) }
      
      # Default palettes if not supplied
      if (is.null(colours_categorical_covariate)) { colours_categorical_covariate <- RColorBrewer::brewer.pal(max(3, length(all_levs)), "Set2") }
      if (is.null(colours_numerical_covariate)) { colours_numerical_covariate <- c("#004643", "#b8b8ff", "#064789") }
      
      # Colour mappinng
      categorical_covariates <- covariates[!covariates %in% numeric_covariates]
      # get all unique categorical levels
      all_levs <- unique(meta_long$Value[meta_long$CovType == "categorical"])
      # assign palette
      merged_colour_map <- setNames( colours_categorical_covariate[seq_along(all_levs)], all_levs )
      
      
      # Plot faceted scatter plot – looping categorical and continuous separately
      {
        ggplot() + 
          # Categorical covariates
          
          geom_point( data = filter(meta_long, CovType == "categorical"),
                      aes(x = PC1, y = PC2, colour = Value),
                      size = 3, alpha = 0.8 ) +
          
          scale_colour_manual( values = merged_colour_map, na.value = "grey50" ) +
          
          # reset the colour scale for continuous
          new_scale_colour() +
          
          # Continuous covariates
          geom_point( data = filter(meta_long, 
                                    CovType == "continuous"),
                      aes(x = PC1, y = PC2, colour = as.numeric(Value)),
                      size = 3, alpha = 0.8 ) +
          
          scale_colour_gradient2(low = colours_numerical_covariate[[1]],
                                 mid = colours_numerical_covariate[[2]],
                                 high = colours_numerical_covariate[[3]],
                                 midpoint = median(as.numeric(filter(meta_long, CovType == "continuous")$Value), na.rm = TRUE),
                                 na.value = "grey50",
                                 name = "Numeric values") +
          
          # Apply facet and common aesthetics
          stat_chull( data = filter(meta_long, CovType == "categorical"),
                      aes(x = PC1, y = PC2, fill = Value, group = Value),
                      alpha = 0.1, geom = "polygon", show.legend = FALSE) +
          facet_wrap(~Covariate, scales = "free") +
          theme_bw(base_size = 14) +
          labs( title = "PCA of Samples",
                subtitle = if(adjusted == TRUE) { paste0("After Covariate Adjustment (", normalisation, ")")} else {paste0("Before Covariate Adjustment (", normalisation, ")")},
                x = "PC1",
                y = "PC2") +
          theme( plot.title = element_text(size =14, face = "bold", hjust = 0.5 ),
                 plot.subtitle = element_text(size =14, face = "italic", hjust = 0.5))
        
      }
      
    } else { message("Transform the matrix: Genes in rows and samples in columns")}
    
  } else { message("matrix or dataframe required as input")}
  
}



# 2. Read all RDS objects from a specified folder
#-------------------------------------------------------------------------------
# Reads all .rds files from a specified directory and combines them into a named 
# list. The supplied directory must contain one or more .rds files.Each RDS file 
# should contain a valid R object. It returns a named list, where each element 
# contains the object read from one RDS file and element name corresponds to the
# original filename.

readRDS_folder <- function(folder){
  
  # Identify all RDS files in the input folder
  files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
  
  # Read each RDS file and store the resulting objects in a list
  out <- lapply(files, function(f) { readRDS(f) })
  
  # Use the corresponding filenames to label the elements of the output list
  names(out) <- basename(files)  
  
  return(out)
}



# 3. Extracting a specified column from multiple data frames
#-------------------------------------------------------------------------------
# Extracts a specified variable from multiple data frames, associates the values
# with gene identifiers, resolves duplicate gene entries, and merges all datasets
# into a single gene-level table.Duplicate entries containing only missing 
# values become one NA value. Duplicate entries with the same value collapse 
# to one value. Duplicate entries with different values retained as all unique values .
# df_list =	Named list of data frames
# column	= Character string specifying the column to extract
# It returns a merged data frame containing genes occurring in only some datasets 
# are retained, with missing values represented by NA.

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



# 4. Helper: single-gene DL meta-analysis
#-------------------------------------------------------------------------------
# Performs DerSimonian–Laird meta-analysis for a single gene using effect 
# estimates and their standard errors from multiple datasets. It supports 
# either fixed-effect or random-effects models.
# beta = Numeric vector of effect estimates
# se = Numeric vector of standard errors
# The two vectors must correspond element-by-element. The function removes
# non-finite effect estimates and SEs before analysis. Returns a named list.

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



# 5. Vectorized wrapper: gene-matrix DL meta-analysis
#-------------------------------------------------------------------------------
# Applies the DerSimonian–Laird meta-analysis independently to every gene in an 
# expression/differential-expression result matrix. It therefore extends the 
# single-gene meta-analysis to genome-wide or multi-gene datasets. Rows of both 
# matrices must correspond to the same genes. Columns must correspond to the 
# same datasets and be in the same order.
# logFC_mat	= Matrix of gene-level effect estimates
# SE_mat	= Matrix of corresponding standard errors
# genes	= Optional gene names
# method	= "random" or "fixed"
# min_se	= Minimum allowed SE
# p_adjust	= "none", "BH" or "bonferroni"
# It returns a dataframe containing gene, k, beta, se, z, p, p_adj, Q, tau2, I2.

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
    out <- meta_dl_single_gene( beta = beta_i,
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



# 6. Export edge files of a network
#-------------------------------------------------------------------------------
# Extracts a selected gene network from a WGCNA TOM matrix, evaluates how network 
# density changes across TOM thresholds, and optionally exports the network in a
# Cytoscape-compatible format. The ggplot component shows the relationship between 
# TOM threshold and number of retained edges. The optional cytoscape component 
# contains the Cytoscape export result.
# geneset = 	Module name or vector of genes; TOM = WGCNA TOM matrix;
# blockwiseModules_network = 	WGCNA module object required for module mode;
# GeneAnnotation = 	Optional mapping dataframe; 
# geneset_type = 	"module" or "genes"
# thresholds =	TOM thresholds to evaluate; 
# export =	Logical; whether to export Cytoscape files
# export_threshold =	TOM threshold used for Cytoscape export
# prefix	 = Output filename prefix
# file_format = 	.txt, .tsv, or .csv

GenesetCytoscapeExport <- function(geneset,
                                   TOM,
                                   blockwiseModules_network = NULL,
                                   GeneAnnotation = NULL,
                                   geneset_type = c("module", "genes"),
                                   thresholds = seq(0, 0.1, by = 0.01),
                                   export = FALSE,
                                   export_threshold = NULL,
                                   prefix = "network",
                                   file_format = c(".txt", ".tsv", ".csv")) {
  
  # Sets input options
  geneset_type <- match.arg(geneset_type)
  file_format <- match.arg(file_format)
  
  # Identifies genes from a WGCNA module or user-defined gene set
  if (geneset_type == "module") {
    
    if (is.null(blockwiseModules_network))
      stop("blockwiseModules_network is required for module analysis.")
    
    inModule <- blockwiseModules_network$colors == geneset
    genes <- colnames(TOM)[inModule]
    
  } else {
    
    genes <- intersect(geneset, rownames(TOM))
    inModule <- rownames(TOM) %in% genes
  }
  
  # Displays message if no adjacency is found for the genes
  if (length(genes) < 2)
    stop("Fewer than 2 genes found in TOM.")
  
  # Extracts the TOM subnetwork corresponding to the selected genes
  modTOM <- TOM[genes, genes, drop = FALSE]
  
  # Cross maps using user-defined annotation table
  if (!is.null(GeneAnnotation)) {
    map <- GeneAnnotation[match(rownames(modTOM), GeneAnnotation$ENSEMBL), ,drop = FALSE]
    keep <- !is.na(map$SYMBOL) & map$SYMBOL != ""
    modTOM <- modTOM[keep, keep, drop = FALSE]
    genes <- make.unique(map$SYMBOL[keep])
    rownames(modTOM) <- genes
    colnames(modTOM) <- genes
  }
  
  # Calculates pairwise TOM values and edge counts across thresholds
  tomValues <- modTOM[upper.tri(modTOM, diag = FALSE)]
  edge_counts <- sapply(thresholds,\(x) sum(tomValues > x, na.rm = TRUE))
  edge_table <- data.frame(Threshold = thresholds,
                           Edges = edge_counts)
  print(edge_table)
  
  # Plots the no.of edges found at respective TOM value threshold
  p <- ggplot(edge_table, aes(Threshold, Edges)) +
    geom_point() +
    geom_line() +
    labs(x = "Threshold",
         y = "Edges at threshold") +
    theme_classic()
  
  # Optional: Export the gene network to Cytoscape when requested 
  cytoscape <- NULL
  if (export) {
    if (is.null(export_threshold))
      stop("export_threshold must be supplied when export = TRUE.")
    cytoscape <- exportNetworkToCytoscape(modTOM,
                                          edgeFile = paste0("CytoscapeInput-edges-", prefix, file_format),
                                          nodeFile = paste0("CytoscapeInput-nodes-", prefix, file_format),
                                          weighted = TRUE,
                                          threshold = export_threshold,
                                          nodeNames = colnames(modTOM),
                                          altNodeNames = colnames(modTOM),
                                          nodeAttr = modTOM,
                                          includeColNames = TRUE)
  }
  
  # Stores the selected genes, TOM network, threshold analysis, plot, and export results
  list(inModule = inModule,
       genes = genes,
       modTOM = modTOM,
       tomValues = tomValues,
       thresholds = thresholds,
       edge_counts = edge_counts,
       edge_table = edge_table,
       ggplot = p,
       cytoscape = cytoscape)
}



# 7. Crossmap
#-------------------------------------------------------------------------------
# Converts gene identifiers between Ensembl, HGNC symbol, and Entrez Gene 
# identifiers using the locally stored gene annotation table. The annotation RDS
# file must exist at the path specified inside the function. If id_col = NULL, 
# the function uses row names as the source identifiers. Returns the same 
# data frame with row names replaced by the mapped identifiers. Returns the same 
# data frame with row names replaced by the mapped identifiers. Unmapped genes 
# receive identifiers of the form: UNMAPPED_originalID
# Duplicate mapped identifiers are made unique using make.unique().
# df = Data frame containing genes
# id_col = Optional column containing source identifiers
# from = Source ID type: "ensembl", "symbol", "entrez"
# to = Target ID type

crossmap <- function(df,
                     id_col = NULL,
                     from = c("ensembl", "symbol", "entrez"),
                     to   = c("ensembl", "symbol", "entrez")) {
  # Loads gene annotation table
  gene_table <- readRDS( "D:\\DRYLAB\\ATOPICDERMATITIS\\Mapping\\GeneAnnotation_ensembl.v86_corrected.rds")
  
  # Sets source and target identifier types
  from <- match.arg(from)
  to   <- match.arg(to)
  
  # Defines corresponding annotation columns
  col_map <- list(ensembl = "ENSEMBL",
                  symbol  = "SYMBOL",
                  entrez  = "ENTREZID")
  from_col <- col_map[[from]]
  to_col   <- col_map[[to]]
  
  # Extracts identifiers from the specified column or row names
  if (!is.null(id_col)) {
    old_ids <- as.character(df[[id_col]])
  } else {
    old_ids <- rownames(df)
  }
  
  # Standardizes identifier types for matching
  gene_table[[from_col]] <- as.character(gene_table[[from_col]])
  gene_table[[to_col]]   <- as.character(gene_table[[to_col]])
  
  # Maps source identifiers to target identifiers
  new_ids <- gene_table[[to_col]][match(old_ids, gene_table[[from_col]])]
  
  # Reports and retains unmatched identifiers
  unmatched <- is.na(new_ids)
  cat("Mapped:", sum(!unmatched), "/", length(new_ids), "\n")
  cat(old_ids[is.na(new_ids)])
  
  # Assigns unique row names while preserving unmatched identifiers
  new_ids[unmatched] <- paste0("UNMAPPED_", old_ids[unmatched])
  
  rownames(df) <- make.unique(new_ids)
  
  return(df)
}



# 8. Pooling recurrent hub genes across CytoHubba ranking methods
#-------------------------------------------------------------------------------
# Identifies genes that repeatedly occur among multiple CytoHubba hub-gene 
# ranking methods within each module/ genesets. Each CytoHubba result data frame must 
# contain a column named: Name. group_list should be a nested list structured 
# approximately as: Module/ genesets
# ├── Method1 → data.frame
# ├── Method2 → data.frame
# ├── Method3 → data.frame
# └── ...
# group_listv = Nested list of CytoHubba ranking data frames
# cutoff = Minimum number of ranking lists in which a gene must occur
# Returns a nested list with the same module/group organization. For each module,
# the output contains: Genes, Freq ; where Freq represents the number of CytoHubba
# ranking lists containing the gene. If no gene reaches the requested cutoff, 
# the function progressively lowers the cutoff until at least one gene is identified.

Get_Common_HUB_Genes_flexible <- function(group_list, cutoff){
  
  # Defines progressively lower cutoffs if the primary cutoff is not met
  altcutoffs <- 1:(cutoff-1)
  
  # Recursive function for processing nested module/ranking lists
  process_node <- function(node, path = "") {
    
    # Leaf node: processes multiple CytoHubba ranking data frames
    if (is.list(node) &&
        length(node) > 0 &&
        all(sapply(node, is.data.frame))) {
      
      cat("\n--- Processing:", path, "---\n")
      
      # Extracts hub-gene names from all ranking methods
      genes <- as.character(unlist(lapply(node, function(df) df$Name)))
      cat("Extracted", length(genes), "genes\n")
      
      # Counts the number of ranking lists in which each gene occurs
      module <- table(genes)
      module <- module[module >= cutoff]
      module <- as.data.frame(module)
      
      # Adaptive cutoff when the predefined threshold is too stringent
      if (nrow(module) == 0) {
        
        message("Cutoff ", cutoff," too high for ", path, " — trying lower cutoffs.")
        
        for (a in altcutoffs) {
          altcutoff <- cutoff - a
          genes_tmp <- unlist(lapply(node, function(df) df$Name))
          module_tmp <- table(genes_tmp)
          module_tmp <- module_tmp[module_tmp >= altcutoff]
          module_tmp <- data.frame(Genes = names(module_tmp),
                                   Freq = as.vector(module_tmp),
                                   stringsAsFactors = FALSE)
          if (nrow(module_tmp) > 0) {
            message("Using cutoff ", altcutoff,
                    " for ", path)
            module <- module_tmp
            break
          }
        }
      }
      
      # Standardizes output column names  
      if (ncol(module) == 2)
        colnames(module) <- c("Genes","Freq")
      
      return(module)
    }
    
    # Nested node: recursively processes each module/ranking group
    if (is.list(node)) {
      
      out <- list()
      
      for (nm in names(node)) {
        
        new_path <- if (path == "")
          nm else paste(path, nm, sep = "/")
        
        out[[nm]] <- process_node(node[[nm]], new_path)
      }
      
      return(out)
    }
    
    # Returns NULL for unsupported input structures
    return(NULL)
  }
  
  # Runs recursion
  return(process_node(group_list))
}



# 9. Hub gene analysis using igraph package
#-------------------------------------------------------------------------------
# Constructs gene interaction networks from Cytoscape edge files, calculates 
# multiple network-centrality measures, combines selected measures into a hub 
# score, identifies top hub genes, and generates network visualizations. The edge 
# file should contain at least two columns representing source and target genes
# and preferably a weight column. The hub score represents the sum of standardized 
# centrality measures, not a single native network-centrality metric.
# modules = characte vector of names of intrested modules/ genesets,
# edge_dir = Directory containing Cytoscape edge files,
# output_dir = Directory for results,(Degree/ WeightedDegree/ Betweenness/
#              Closeness/ Eigenvector/ PowerCentrality/PageRank)
# methods = Centrality measures used for hub scoring,
# hubcolor = hexcode to form gradient with color for hub gene at an extreme   
# nodecolor = hexcode to form gradient with color for node gene at an extreme   
# top_n = Number of top hub genes
# plot_type = "igraph" or "ggraph",
# export_type = "png" or "pdf"
# Returns a list indexed by module. Each module contains:centrality, top_genes, graph
analyze_hub_modules <- function(modules,
                                edge_dir,
                                output_dir,
                                methods,
                                hubcolor,
                                nodecolor,
                                top_n,
                                plot_type = c("igraph", "ggraph"),
                                export_type = c("png", "pdf")) {
  suppressPackageStartupMessages({
    library(igraph)
    library(ggraph)
    library(tidygraph)
    library(ggrepel)
    library(RColorBrewer)
    library(scales)
    library(dplyr)
    library(viridis)
    library(ggplot2)
  })
  
  # Sets the plot type and output file format
  plot_type <- match.arg(plot_type)
  export_type <- match.arg(export_type)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  all_results <- list()
  g_tbl_nodes <- list()
  g_tbl_edges <- list()
  
  # Defines a continuous colour gradient for hub score
  col_fun <- colorRampPalette(c(nodecolor, hubcolor)) # pink to red gradient col_fun <- colorRampPalette(c("#4575B4", "#D73027", "#4575B4")) # red to blue gradient
  
  for (mod in modules) {
    message("Processing module: ", mod)
    
    # Loads the exported cytoscape edge file
    edge_path <- file.path(edge_dir, paste0("CytoscapeInput-edges-", mod, ".tsv"))
    if (!file.exists(edge_path)) {  warning("Edge file not found for module: ", mod)
      next}
    edges <- read_delim(edge_path, col_names = TRUE)
    
    # Check: graph construction needs more than two nodes/genes
    if(ncol(edges) < 2)
      stop("Invalid edge file for ", mod)
    
    # Contructs the undirected gene interaction network
    # Hub calculation
    {
      g <- graph_from_data_frame(edges, directed = FALSE)
      
      # Restricts the network only to the strongest weighted interactions
      if (vcount(g) > 2000) {
        message("Graph too large (", vcount(g), " nodes). Sampling top edges by weight...")
        edges <- as.data.frame(get.edgelist(g))
        edges$weight <- E(g)$weight
        edges <- edges[order(-edges$weight), ][1:5000, ]  # keep top 5000 strongest edges
        g <- graph_from_data_frame(edges, directed = FALSE)
      }
      
      # Calculates network centrality measures
      deg  <- degree(g)
      wdeg <- strength(g, weights = E(g)$weight)
      bet  <- tryCatch(betweenness(g, weights = 1/E(g)$weight, normalized = TRUE), error = function(e) rep(0, vcount(g)))
      clo  <- tryCatch(closeness(g, weights = 1/E(g)$weight, normalized = TRUE),  error = function(e) rep(0, vcount(g)))
      eig  <- tryCatch(eigen_centrality(g, weights = E(g)$weight)$vector,         error = function(e) rep(0, vcount(g)))
      powcen <- tryCatch(power_centrality(g, rescale = FALSE),                    error = function(e) rep(0, vcount(g)))
      pr   <- tryCatch(page_rank(g, weights = E(g)$weight)$vector,                error = function(e) rep(0, vcount(g)))
      
      centrality_table <- data.frame(Gene = names(deg),
                                     Degree = deg,
                                     WeightedDegree = wdeg,
                                     Betweenness = bet,
                                     Closeness = clo,
                                     Eigenvector = eig,
                                     PowerCentrality = powcen,
                                     PageRank = pr,
                                     stringsAsFactors = FALSE)
      
      # Standardizes centrality method names supplied by the user
      methods[methods %in% c("degree", "deg", "Degree")] <- "Degree"
      methods[methods %in% c("Betweeness", "Betweenness", "betweenness", "between", "bet", "Bet")] <- "Betweenness"
      methods[methods %in% c("powcen", "power", "pc", "PowerCentrality", "powerCentrality")] <- "PowerCentrality"
      methods[methods %in% c("pagerank", "pr", "PageRank")] <- "PageRank"
      methods[methods %in% c("WeightedDegree", "wdeg", "weightedDegree")] <- "WeightedDegree"
      methods[methods %in% c("Closeness", "clo", "closeness")] <- "Closeness"
      
      # Check: method specified not in centrality measures
      missing_methods <- methods[!methods %in% colnames(centrality_table)]
      if(length(missing_methods) > 0){
        stop("These centrality methods are invalid or misnamed: ",
             paste(missing_methods, collapse = ", "),
             "\nValid options are: ",
             paste(colnames(centrality_table)[-1], collapse=", "))
      }
      
      # Calculates the combined hub score from selected centrality measures
      centrality_table$HubScore <- rowSums(scale(centrality_table[ ,methods, drop = FALSE]), na.rm = TRUE)
      
      # Rank genes according to the hub score
      centrality_table <- centrality_table[order(-centrality_table$HubScore), ]
      topHUB_genes <- head(centrality_table$Gene, top_n)
      
      # Saves the centrality measures and hub scores
      out_csv <- file.path(output_dir, paste0("HubGeneCentrality_", mod, "_", paste(methods, collapse = "_"), ".csv"))
      write.csv(centrality_table, out_csv, row.names = FALSE)
      }
    
    # Maps the hub scores to node colours 
    hubscore_scaled <- rescale(centrality_table$HubScore, to = c(0,1))
    color_values <- col_fun(100)[as.numeric(cut(hubscore_scaled, 100))]
    
    # Prepares the network to tidygraph format
    HubScore <- centrality_table$HubScore[match( V(g)$name, centrality_table$Gene )]
    g_tbl <- as_tbl_graph(g) %>%
      mutate(Module = mod,
             isHub = name %in% topHUB_genes,
             HubScore = HubScore,
             node_size = rescale(HubScore, to = c(3, 15)), 
             node_color = color_values[match(name, centrality_table$Gene)])
    
    # Stores node and edge data ggraph 
    g_tbl_nodes[[mod]] <- g_tbl %>% activate(nodes) %>% as_tibble() %>% mutate(Module = mod)
    g_tbl_edges[[mod]] <- g_tbl %>% activate(edges) %>% as_tibble() %>% mutate(Module = mod)
    
    # igraph plot
    if (plot_type == "igraph" || NULL) {
      
      message("Generating individual igraph plots for each module...")
      set.seed(123)
      
      # output file for the network graph (pdf/png)
      outfile <- file.path(output_dir, paste0("HubGeneCentrality_", mod, "_", paste(methods, collapse = "_"), "igraph.", export_type))
      if (export_type == "pdf") {
        pdf(outfile, width = 8, height = 8, useDingbats = FALSE)
      } else {
        png(outfile, width = 2500, height = 2500, res = 350, bg = "transparent", antialias = "cleartype")
      }
      
      layout <- layout_with_fr(g, niter = 2000)
      layout <- layout * 2 
      
      hub_colors <- col_fun(100)                # Assigns a hub score
      
      node_colors <- sapply(V(g)$name, function(name) { 
        if (name %in% topHUB_genes) { color_values[match(name, centrality_table$Gene)]
        } else {"#0047AB"}
      })
      
      label_colors <- ifelse(V(g)$name %in% topHUB_genes, "black", "black")  
      
      vertex_shapes <- ifelse(V(g)$name %in% topHUB_genes, "vrectangle", "circle")
      
      # Plot
      plot( g,
            layout = layout,
            vertex.color = node_colors,
            vertex.shape = vertex_shapes,
            vertex.size = ifelse(V(g)$name %in% topHUB_genes, 16, 3),
            vertex.label = ifelse(V(g)$name %in% centrality_table$Gene, V(g)$name, NA),
            vertex.size2 = ifelse(V(g)$name %in% topHUB_genes, 5,3),
            vertex.label.cex = ifelse(V(g)$name %in% topHUB_genes, 0.6, 0.5),
            vertex.label.color = label_colors,
            vertex.frame.color = "white",
            edge.color = alpha("black", 0.3),
            main = paste("Hub Network -", mod),
            asp = 1,
            bg = "white"
      )
      par(xpd = TRUE)                            
      legend( "top",                           
              legend = c("Other Genes","Top 10 Hub GeneS"),
              pch = c(21, 22),
              col = "black",
              pt.bg = c("#0047AB", "white"),
              pt.cex = 1.5,
              bty = "n",
              cex = 0.7,
              text.font = 2,
              ncol = 2)
      # Hub score gradient legend
      legend_x <- -0.026           
      symbol_w <- strwidth("Top 10 Hub Genes  ", cex = 0.7) * 0.6      
      symbol_h <- strheight("Top 10 Hub Genes  ", cex = 0.7) * 1.8
      grad_left   <- legend_x - symbol_w * 0.3
      grad_right  <- legend_x + symbol_w * 0.3
      grad_top    <- legend_y + symbol_h * 0.35
      grad_bottom <- legend_y - symbol_h * 0.35
      gradient_colors <- col_fun(100)
      rect_x <- seq(grad_left, grad_right, length.out = length(gradient_colors) + 1)
      for (i in seq_along(gradient_colors)) {
        rect(rect_x[i], grad_bottom, rect_x[i + 1], grad_top,
             col = gradient_colors[i], border = NA)
      }
      par(xpd = FALSE)
      dev.off()
    }
    
    # Stores centrality results and network object the module
    all_results[[mod]] <- list(centrality = centrality_table,
                               top_genes = topHUB_genes,
                               graph = g)
    
  }
  
  # Individual ggraph plots per module
  if (plot_type == "ggraph") {
    message("Generating individual ggraph plots for each module...")
    for (mod in names(g_tbl_nodes)) {
      message("Plotting module: ", mod)
      
      g_tbl_nodes_mod <- g_tbl_nodes[[mod]]
      g_tbl_edges_mod <- g_tbl_edges[[mod]]
      
      # Adjusts hub-node size according to gene-label length
      g_tbl_nodes_mod <- g_tbl_nodes[[mod]] %>%
        mutate(label_length = nchar(name),
               node_size = ifelse( isHub, pmin(35, 8 + label_length * 1.5), 4 ),
               hub_status = ifelse(isHub, "Hub gene", "Non-hub gene"))
      g_tbl_edges_mod <- g_tbl_edges[[mod]]
      g_tbl_mod <- tbl_graph( nodes = g_tbl_nodes_mod,
                              edges = g_tbl_edges_mod,
                              directed = FALSE )
      V(g_tbl_mod)$NodeType <- ifelse(V(g_tbl_mod)$isHub, "Hub gene", "Non-hub gene")
      
      # Plot
      {
        set.seed(123)
        g_plot <- ggraph(g_tbl_mod, layout = "fr", niter = 2000) +
          geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.4, show.legend = TRUE) +
          scale_edge_width(name = "Edge Weight", range = c(0.3, 2.0)) +
          geom_node_point( aes(fill = HubScore, size = NodeType),
                           shape = 21,
                           color = "gray24",
                           stroke = 0.4,
                           show.legend = TRUE) +
          scale_size_manual( name = "Node Type",
                             values = c("Hub gene" = 15, "Non-hub gene" = 2),
                             guide = guide_legend( title.position = "top", 
                                                   title.hjust = 0.5, 
                                                   override.aes = list(shape = 21)) )+
          # Hub gene label centered in the node
          geom_node_text(data = function(x) dplyr::filter(x, isHub),
                         aes(label = name),
                         color = "black",
                         fontface = "bold",
                         size = 4,
                         vjust = 0.5, hjust = 0.5,
                         lineheight = 0.9, show.legend = FALSE) +
          
          # Repel gene label centered in the node
          geom_node_text(data = function(x) dplyr::filter(x, !isHub),
                         aes(label = name),
                         repel = TRUE,
                         color = "black", size = 2.7,
                         fontface = "bold",
                         box.padding = 0.2, point.padding = 0.1,
                         segment.color = "gray40",
                         segment.size = 0.15,
                         max.overlaps = 15, show.legend = TRUE) +
          
          # Hub score color legend
          scale_fill_gradientn( colours = col_fun(100),
                                name = "Hub Score",
                                limits = c(min(V(g_tbl_mod)$HubScore, na.rm = TRUE),
                                           max(V(g_tbl_mod)$HubScore, na.rm = TRUE)),
                                guide = guide_colorbar(barwidth = 5, barheight = 0.5,
                                                       title.position = "top",
                                                       title.hjust = 0.5)) +
          theme_void() +
          theme( legend.position = "bottom",
                 legend.box = "horizontal",
                 legend.title = element_text(face = "bold", size = 7),
                 legend.text = element_text(size = 6),
                 plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                 plot.background = element_rect(fill = "white", color = NA)) +
          ggtitle(paste("Hub Gene Network -", mod))
      }
      
      # Saves each plot individually
      outfile <- file.path(output_dir, paste0("HubGeneCentrality_", mod, "_", paste(methods, collapse = "_"), "ggraph.", export_type))
      ggsave(outfile, plot = g_plot, width = 8, height = 6, dpi = 350, bg = "transparent")
    }
  }
  
  return(all_results)
  
}



# 10. Hub gene analysis using intramodularConnectivity()
#-------------------------------------------------------------------------------
# Identifies and ranks WGCNA hub genes using intramodular connectivity together 
# with module membership. It can analyse either existing WGCNA modules or 
# user-defined gene sets. Returns a named list containing one data frame per module.
# Genes are ranked primarily by kWithin and then by module membership.
# modules = character vector of names of interested modules 
# module_membership_measure = dataframe of modules membership
# adjacency = adjacency matrix
# wgcna_colors = character vector generated by blockwise module having genes and
#                respective module color
# GeneSymbols = dataframe of gene mapping having ensemble-ID, entrezID, and hgnc_symbol
# gene_sets = a simple list of character vectors containing gene ensemble-ID 
#             (optional if the name of the module is not present wgcna_colors, else NULL)
HUB_intramodularConnectivity <- function(modules, 
                                         module_membership_measure, 
                                         adjacency, 
                                         wgcna_colors,
                                         GeneSymbols,
                                         gene_sets) {
  # Specifies output file format
  file_type <- readline(prompt = "Enter desired output format (.csv or .xlsx): ")
  file_type <- tolower(trimws(file_type))
  
  # Calculates intramodular connectivity for all genes
  cat("\n Calculating intramodularConnectivity")
  IMConn <- intramodularConnectivity(adjacency, wgcna_colors)
  IMConn$Genes <- rownames(IMConn)
  View(IMConn)
  
  hub_list <- list()
  
  for (module in modules) {
    
    # Module specific calculation
    cat(paste("\n Calculating intramodularConnectivity for module", module))
    
    # Idntifies genes belonging to a module or separate geneset
    if(module %in% wgcna_colors){
      cat(paste("\n Fetching intramodularConnectivity for module", module))
      module_hubs <- IMConn[wgcna_colors == module, ]
      module_hubs$module.membership <- module_membership_measure[wgcna_colors == module, paste0("ME", module)]
    }
    else {
      cat(paste("\n WGCNA color not found for module", module, "\n Fetching best match based on highest module membership"))
      if (is.null(gene_sets)) {
        stop("\nERROR: item is not a module and no gene_sets provided.\n")
      }
      if (!module %in% names(gene_sets)) {
        stop("\nERROR: item is neither module nor a gene_set name.\n")
      }
      genes <- gene_sets[[module]]
      genes_present <- intersect(genes, rownames(IMConn))
      length(genes_present)
      if (length(genes_present) == 0) {
        warning(paste("\nNo genes from", module, "found in IntramodularConnectivity. Skipping."))
        next
      }
      
      module_hubs <- IMConn[genes_present, ]
      
      # Extracts only numeric ME columns
      kme_sub <- module_membership_measure[genes_present, ]
      kme_sub <- kme_sub[, sapply(kme_sub, is.numeric), drop = FALSE]
      
      # Computes max module membership per gene
      max_kME <- apply(kme_sub, 1, max, na.rm = TRUE)
      
      # Identifies module with highest membership
      best_module <- colnames(kme_sub)[max.col(kme_sub, ties.method = "first")]
      
      module_hubs$module.membership <- max_kME
      module_hubs$module.name <- best_module
    }
    
    # Stores number of rows before further joins
    rows <- nrow(module_hubs)
    
    # Map gene identifiers to gene symbols and rank hub genes
    cat(paste("\n Consistent Mapping for", module))
    module_hubs_ranked <- module_hubs %>%
      left_join(GeneSymbols, by = c("Genes" = "ENSEMBL")) %>%
      mutate(module.name = as.character(module.name)) %>%
      select(Genes, SYMBOL,
             module.membership, module.name,
             kWithin, kDiff, kTotal, kOut)
    
    # Ranks genes by intramodular connectivity and module membership
    module_hubs_ranked <- module_hubs_ranked[1:rows, ] %>%
      arrange(desc(kWithin), desc(module.membership)) %>%
      distinct(SYMBOL, .keep_all = TRUE)
    
    # Reassign row names as gene symbols
    rownames(module_hubs_ranked) <- module_hubs_ranked$SYMBOL
    
    # Export hub-gene results in the selected format
    if (file_type == ".csv") {
      write.csv(module_hubs_ranked, paste0(module, "_HubGenes.csv"), row.names = TRUE)
    } else if (file_type == ".xlsx") 
    {
      write.xlsx(module_hubs_ranked, paste0(module, "_HubGenes.xlsx"), rowNames = TRUE)
    } else {
      warning("Invalid file format specified. Skipping file export.")
    }
    
    # Stores ranked hub-gene results for each module or gene set
    hub_list[[module]] <- module_hubs_ranked
  }
  
  return(hub_list)
}



# 11. Plot heat map for skin conditions
#-------------------------------------------------------------------------------
# Generates a heat map of selected genes across samples and annotates samples 
# according to skin condition (Healthy, Non-lesional, Lesional). Expression values 
# are Z-score standardized gene-wise across samples, then restricted to [-3, 3] 
# for visualization.
# FullExp = 	Expression matrix
# genes =	Character vector of genes
# Metadata =	Sample metadata
# ann_colors =	Optional annotation colours
# clustering_method =	Hierarchical clustering method
# clustering_distance_rows =	Distance measure for genes
# clustering_distance_cols =	Distance measure for samples
# Plot_title =	Heat-map title
# colorRampPalette =	Optional colour palette
plotheatmap <- function( FullExp,
                         genes, 
                         Metadata, 
                         ann_colors = NULL, 
                         clustering_method = "ward.D2",
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "correlation",
                         Plot_title,
                         colorRampPalette = NULL){
  
  # Z-score normalize expression values across samples for each gene
  filZ <- t(scale(t(FullExp[genes, ])))
  # Limits the normalised expression values to the range [-3, 3] to reduce
  # the influence of extreme values on heat-map visualization
  filZ_clip <- pmin(pmax(filZ, -3), 3)
  range(filZ_clip)
  
  # Generates sample annotation based on disease condition
  fil_htmp_ann  <- data.frame(Condition = Metadata$Condition)
  rownames(fil_htmp_ann ) <- colnames(filZ_clip)
  fil_htmp_ann$Condition <- factor( fil_htmp_ann $Condition, levels = c("Healthy", "Non-lesional", "Lesional"))
  
  # Arrange samples according to disease condition
  ord <- order(fil_htmp_ann$Condition)
  filZ_clip <- filZ_clip[, ord]
  fil_htmp_ann <- fil_htmp_ann [ord, , drop = FALSE]
  
  # Defines default annotation colors for disease conditions
  if(is.null(ann_colors)){ann_colors <- list( Condition = c( Healthy = "#4F6367", `Non-lesional` = "#84BD00", Lesional = "#0077C8"))}
  
  # Defines the default diverging color scale for standardized expression
  if(is.null(colorRampPalette)){
    heat_colors <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)} else {
      heat_colors <- colorRampPalette}
  
  # Applies default clustering parameters when not explicitly specified
  if(is.null(clustering_method)){clustering_method = "ward.D2"}
  if(is.null(clustering_distance_rows)){clustering_distance_rows = "correlation"}
  if(is.null(clustering_distance_cols)){clustering_distance_cols = "correlation"}
  
  # Generates heat map with gene clustering and condition annotations
  fil_htmp <- pheatmap(
    filZ_clip,
    annotation_col = fil_htmp_ann,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    clustering_method = clustering_method,
    clustering_distance_rows = clustering_distance_rows,
    clustering_distance_cols = clustering_distance_cols,
    color = heat_colors,
    fontsize = 16,
    main = Plot_title )
  
  # Renders the heat map in a new graphics device
  dev.new()
  grid::grid.draw(fil_htmp$gtable)
}



# 12. PCA plot
#-------------------------------------------------------------------------------
# Performs PCA on a selected gene-expression matrix and generates a two-dimensional 
# PCA plot showing sample distribution according to skin condition.
# pcaplot ggplot PCA visualization; var_exp = variance explained by each PC;
# pcascores = sample PCA coordinates; pca = complete prcomp object
pca_val <- function(Expression, Metadata, colour_guide, PCs = c("PC1", "PC2"), title){
  
  # Performs principal component analysis on samples
  # Expression matrix is transposed to obtain samples × genes format
  pca <- prcomp(t(Expression))
  
  # Extracts PCA scores for individual samples
  pcascores <- as.data.frame(pca$x)
  pcascores$Sample <- rownames(pcascores)
  
  # Maps sample identifiers to metadata and assigns disease condition
  pcascores$Condition <- Metadata$Condition[match(pcascores$Sample, Metadata$geo_accession)]
  
  # Calculates the proportion of total variance explained by each principal
  # component
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  
  # Definig the x- axis and yacis label representing the percentage % of variance
  xlab <- paste0(PCs[1], " (", round(100 * var_exp[1], 1), "%)")
  ylab <- paste0(PCs[2], " (", round(100 * var_exp[2], 1), "%)")
  
  # Generates a 2D scatter plot of PC1 and PC2
  # Defines colors to use for skin conditions
  if(is.null(colour_guide)){
    cols <- c( "Healthy" = "#4f6367", "Lesional" = "#0077c8", "Non-lesional" = "#84bd00" )
  } else{cols = colour_guide }
  
  pcaplot <- ggplot(pcascores, aes(x = .data[[PCs[1]]],
                                   y = .data[[PCs[2]]],
                                   colour = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_chull( aes(fill = Condition, group = Condition),
                alpha = 0.15,
                geom = "polygon",
                colour = "black",      
                linewidth = 0.5,
                show.legend = FALSE) +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw(base_size = 14) +
    theme( plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
           axis.title.x    = element_text(size = 16, face = "bold"),
           axis.title.y    = element_text(size = 16, face = "bold"),
           axis.text       = element_text(size = 14),
           legend.position = "top",
           legend.title    = element_blank(),
           legend.text    = element_text(size = 14, face = "bold"),
           panel.border    = element_rect(linewidth = 1),
           plot.margin     = margin(10, 10, 10, 10)
    )
  
  pcaplot
  
  return( list(pcaplot = pcaplot, var_exp = var_exp, pcascores = pcascores, pca = pca ))
  
}



# 13. UMAP plot
#-------------------------------------------------------------------------------
# Performs uniform manifold approximation projection dimensionality reduction
# on a gene-expression matrix and visualizes samples according to skin condition.
# The expression matrix is transposed to samples × genes and UMAP is performed
#using n_neighbors = 15, min_dist = 0.1, metric = "euclidean". The function
# returns UMAP_data containing UMAP coordinates and sample metadata and plot.
plot_umap <- function(expression_data, colour, Metadata, title){
  
  # Transposes expression matrix to obtain samples × genes format
  X <- t(expression_data)
  
  # Matches sample identifiers to metadata and assigns disease condition
  Y <- Metadata$Condition[match(rownames(X), Metadata$geo_accession)] 
  
  # Sets seed to ensure reproducibility of the UMAP
  set.seed(123)
  
  # Generates UMAP  using Euclidean distance
  um <- umap(X, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
  umap_df <- data.frame(UMAP1 = um[, 1],
                        UMAP2 = um[, 2],
                        Sample = rownames(um),
                        Condition = Y)
  
  # UMAP scatter plot colored by disease condition
  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Condition)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = colour) +
    theme_bw(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size =14, face = "bold"),
      legend.position = "top",
      axis.title.x    = element_text(size = 16, face = "bold"),
      axis.title.y    = element_text(size = 16, face = "bold"),
      axis.text       = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
    labs(title = title)
  
  p_umap
  
  return(list(UMAP_data = umap_df,
              Plot = p_umap ))
}



# 14. Operator
#-------------------------------------------------------------------------------
# A helper operator that defines a null-coalescing operator to return a fallback value when the primary 
# value is NULL or empty. a %||% b means: use "a" if available; otherwise use "b"
`%||%` <- function(a,b) if(!is.null(a) && length(a)>0) a else b



# 15. JSON
#-------------------------------------------------------------------------------
# Safely extracts and parses JSON content from different types of HTTP response 
# objects used by R web/API requests.  prevents malformed API responses from 
# directly terminating downstream database queries.
# resp = HTTP response object generated by packages such as httr, httr2, or 
#        compatible clients.
# Returns parsed R object/data structure when valid JSON is available or NULL 
# when the response cannot be extracted or parse
safe_json <- function(resp) {
  
  # Extracts and parses raw response content (works for curl/httr response objects)
  if (!is.null(resp$content)) {
    txt <- rawToChar(resp$content)
    if (nzchar(txt)) {
      out <- tryCatch(jsonlite::fromJSON(txt, flatten = TRUE),
                      error = function(e) NULL)
      return(out)
    }
  }
  
  # Handles httr2 response objects
  if ("httr2_response" %in% class(resp)) {
    out <- tryCatch(httr2::resp_body_string(resp), error = function(e) NULL)
    if (!is.null(out)) {
      return(tryCatch(jsonlite::fromJSON(out, flatten = TRUE),
                      error = function(e) NULL))
    }
  }
  
  # Handles generic response objects containing a body field.
  if ("response" %in% class(resp) && !is.null(resp$body)) {
    txt <- rawToChar(resp$body)
    return(tryCatch(jsonlite::fromJSON(txt, flatten = TRUE),
                    error = function(e) NULL))
  }
  
  return(NULL)
}



# 16. safe_query
#-------------------------------------------------------------------------------
# A helper that executes a database query for a single gene while preventing an error for one 
# gene from stopping the complete drug–gene interaction analysis.
# fn =	Database-specific query function
# gene =	Gene symbol to query
# Returns successful query result returned by fn while the error is printed to 
# the console rather than terminating the analysis.
safe_query <- function(fn, gene) {
  tryCatch(
    fn(gene),
    error = function(e) {
      message("[ERROR] ", gene, " in ", deparse(substitute(fn)), ": ", conditionMessage(e))
      tibble()   
    }
  )
}



# 17. build_db_list
#-------------------------------------------------------------------------------
# Applies a database-specific query function to every gene in a supplied gene 
# list and organizes the results into a named list.
# fn =	Database-specific query function
# genes =	Character vector of gene symbols
build_db_list <- function(fn, genes) {
  results <- setNames(
    lapply(genes, function(g) safe_query(fn, g)),
    genes
  )
  return(results)
}



# 18. Drug Central
#-------------------------------------------------------------------------------
# Retrieves drug–gene interaction records for supplied genes from a locally 
# downloaded DrugCentral interaction dataset. The file must contain a column with 
# gene_claim_name. Returns a data frame containing DrugCentral records whose 
# gene_claim_name matches one of the supplied genes.
query_drugcentral <- function(genes) {
  # Locally stored drug–gene interaction file downloaded from DrugCentral 
  # interaction dataset.
  file = "D:\\DRYLAB\\ATOPICDERMATITIS\\Durg Gene-Interaction\\interactions.tsv"
  df <- readr::read_tsv(file, col_types = cols(.default = "c"))
  df <- df[, colSums(!is.na(df)) > 0]
  df
  
  # Specific column to match with the genes in query
  gene_col <- "gene_claim_name"
  
  # Gene symbols are matched case-insensitively to the gene_claim_name field 
  # and columns containing no observed data are removed prior to filtering.
  genes_up <- toupper(genes)
  df %>% filter(toupper(.data[[gene_col]]) %in% genes_up)
}


# 19. ChEMBL 
#-------------------------------------------------------------------------------
# Query ChEMBL for drug–target activity records associated with each gene 
# Returns a tibble containing fields gene, source_db, drug, interaction_types,
# evidence_sources, pubmed_ids, source_url, notes.
# Internet/API access is required.
query_chembl_single <- function(gene) {
  
  # Searches ChEMBL for target records associated with the input gene.
  resp <- GET("https://www.ebi.ac.uk/chembl/api/data/target/search",
              query = list(q = gene, format = "json"))
  js <- safe_json(resp)
  
  # Returns NULL if no valid target records are identified.
  if (is.null(js) || is.null(js$targets) || !is.data.frame(js$targets) || nrow(js$targets) == 0)
    return(NULL)
  
  # Extracts unique ChEMBL target identifiers.
  ids <- unique(js$targets$target_chembl_id)
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return(NULL)
  
  # Extracts individual fields while handling missing values.
  sx <- function(row, field) {
    if (is.null(row)) return("")
    if (is.data.frame(row) && field %in% names(row)) {
      v <- row[[field]][1]
      if (is.null(v) || length(v) == 0 || is.na(v)) return("")
      return(as.character(v))
    }
    if (is.list(row) && field %in% names(row)) {
      v <- row[[field]]
      if (is.null(v) || length(v) == 0 || is.na(v)) return("")
      return(as.character(v))
    }
    return("")
  }
  
  # Retrieves records for each identified ChEMBL target.
  bind_rows(lapply(ids, function(tid) {
    resp2 <- GET("https://www.ebi.ac.uk/chembl/api/data/activity",
                 query = list(target_chembl_id = tid, limit = 200, format = "json"))
    js2 <- safe_json(resp2)
    
    # Skips targets with unavailable or empty activity records.
    if (is.null(js2) || is.null(js2$activities) || !is.data.frame(js2$activities))
      return(NULL)
    acts <- js2$activities
    if (nrow(acts) == 0) return(NULL)
    
    # Extracts and standardizes drug, activity, evidence, and source information.
    bind_rows(lapply(seq_len(nrow(acts)), function(i) {
      a <- acts[i, ]
      tibble(gene  = gene,
             source_db = "ChEMBL",
             drug  = {
               p1 <- sx(a, "molecule_pref_name")
               p2 <- sx(a, "molecule_chembl_id")
               if (nzchar(p1)) p1 else p2
             },
             interaction_types = sx(a, "standard_type"),
             evidence_sources  = sx(a, "assay_type"),
             pubmed_ids        = sx(a, "src_id"),
             source_url = paste0("https://www.ebi.ac.uk/chembl/compound_report_card/",
                                 sx(a, "molecule_chembl_id")),
             notes = paste("value:", sx(a, "standard_value")))
    }))
  }))
}

# Applies query_chembl_single() to multiple genes and combines the results into
# one data frame. Returns A combined tibble containing ChEMBL drug–gene 
# interaction/activity records.
query_chembl <- function(genes) {
  bind_rows(lapply(genes, query_chembl_single))
}


# 20. OpenTargets
#-------------------------------------------------------------------------------
# Retrieves drugs and clinical candidates associated with a human gene from the
# Open Targets Platform. Gene symbol must be recognized by org.Hs.eg.db. 
# Internet/API access is required.
query_opentargets <- function(gene_symbol) {
  
  # Maps the gene symbol to its Ensembl gene identifier.
  ensg <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                keys = gene_symbol,
                                column = "ENSEMBL",
                                keytype = "SYMBOL",
                                multiVals = "first")
  
  # Returns NULL if no Ensembl identifier is available
  if (is.na(ensg)) return(NULL)
  
  # The GraphQL query to retrieve drugs and clinical candidates 
  # associated with the specified Open Targets gene entry.
  gql <- 'query KD($ensemblId: String!) {
             target(ensemblId: $ensemblId) {
               approvedSymbol
               drugAndClinicalCandidates {
                count
                rows {
                 id
                 maxClinicalStage
                 drug {
                  id
                  name
                  drugType }
                             }
                               }
                                 }
                                   }
                                     '
  # Constructs the GraphQL request body using the Ensembl gene identifier.
  body <- list( query = gql,
                variables = list(ensemblId = ensg))
  
  # Submits the query to the Open Targets GraphQL API
  resp <- httr::POST("https://api.platform.opentargets.org/api/v4/graphql",
                     body = body,
                     encode = "json")
  js <- httr::content(resp)
  
  # Stops execution if the API returns a query or server error.
  if (!is.null(js$errors)) {
    stop(paste(vapply(js$errors, `[[`, "", "message"),
               collapse = "\n"))
  }
  
  # Extracts drug and clinical-candidate records associated with the target.
  rows <- js$data$target$drugAndClinicalCandidates$rows
  
  # Returns NULL when no associated drug records are available.
  if (is.null(rows) || length(rows) == 0)
    return(NULL)
  
  # Standardize the retrieved drug and clinical-stage information.
  dplyr::bind_rows( lapply(rows, function(r) {
    tibble::tibble(gene = gene_symbol,
                   ensembl = ensg,
                   drug_id = r$drug$id,
                   drug_name = r$drug$name,
                   drug_type = r$drug$drugType,
                   clinical_stage = r$maxClinicalStage) }))
  
}


# 21. PharmaGKB
#-------------------------------------------------------------------------------
# Retrieves clinical drug–gene annotations from PharmGKB for a supplied set of 
# gene symbols. The function requires access to PharmGKB downloadable data. 
# If the required files are absent, it downloads: genes.zip and 
# clinicalAnnotations.zip files. Returns a tibble containing gene,drug,
# evidenceLevel, phenotypeCategory, pmid_count, annotation_id, url.
# Internet/API access is required.
query_pharmGKB <- function(gene_symbols, dir = "pharmgkb_data") {
  
  # Creates a local directory for PharmGKB data if it does not exist.
  if (!dir.exists(dir)) dir.create(dir)
  
  # The paths for the required data files retrieved from PharmGKB
  # and stored locally
  genes_file <- file.path(dir, "genes.tsv")
  ca_file    <- file.path(dir, "clinical_annotations.tsv")
  
  # Downloads and extract PharmGKB data files when they are not 
  # already available locally.
  if (!file.exists(genes_file) || !file.exists(ca_file)) {
    message("Downloading PharmGKB files...")
    
    files <- c("genes.zip", "clinicalAnnotations.zip")
    base_url <- "https://api.pharmgkb.org/v1/download/file/data/"
    
    for (f in files) {
      url <- paste0(base_url, f)
      dest <- file.path(dir, f)
      
      message("Downloading: ", f)
      tryCatch(
        download.file(url, destfile = dest, mode = "wb"),
        error = function(e) stop("Failed to download: ", f)
      )
      unzip(dest, exdir = dir)
    }
    
    message("Download & extraction complete.")
  }
  
  # Loads PharmGKB gene annotations and retains gene symbols and corresponding 
  # PharmGKB identifiers
  genes <- vroom::vroom( genes_file,
                         delim = "\t",
                         col_types = vroom::cols(.default = "c"),
                         trim_ws = TRUE,
                         progress = FALSE) %>%
    dplyr::select(geneSymbol = Symbol,
                  pharmgkbId = `PharmGKB Accession Id`)
  
  # Restricts the analysis to the predefined genes of interest.
  gene_map <- genes %>%
    dplyr::filter(geneSymbol %in% gene_symbols)
  
  # Returns an empty result if none of the target genes are represented in the  
  # PharmGKB gene annotation file.
  if (nrow(gene_map) == 0) {
    warning("No supplied genes found in PharmGKB genes.tsv")
    return(tibble::tibble())
  }
  
  # Load PharmGKB clinical annotation records.
  ca <- vroom::vroom(ca_file,
                     delim = "\t",
                     col_types = vroom::cols(.default = "c"),
                     trim_ws = TRUE,
                     progress = FALSE)
  
  # Extracts and standardizes drug–gene clinical annotations 
  results <- ca %>%
    dplyr::filter(Gene %in% gene_map$geneSymbol) %>%
    tidyr::separate_rows(`Drug(s)`, sep = ";") %>%
    dplyr::mutate(`Drug(s)` = stringr::str_trim(`Drug(s)`)) %>%
    dplyr::transmute(gene = Gene,
                     drug = `Drug(s)`,
                     evidenceLevel = `Level of Evidence`,
                     phenotypeCategory = `Phenotype Category`,
                     pmid_count = `PMID Count`,
                     annotation_id = `Clinical Annotation ID`,
                     url = URL) %>%
    dplyr::distinct()
  
  return(results)
}


# 22. KEGG DBGET
#-------------------------------------------------------------------------------
# Retrieves drug associations for human genes from the KEGG database using the 
# KEGG REST API. Gene symbol must be recognized by org.Hs.eg.db.  Internet/API 
# access is required.
query_kegg <- function(gene){
  
  # Applies the function iteratively when multiple genes are provided.
  if (length(gene) > 1) {
    return(purrr::map_df(gene, query_kegg))
  }
  
  # Maps the gene symbol to its Entrez Gene identifier.
  entrez <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db,
                                                   keys = gene,
                                                   column = "ENTREZID",
                                                   keytype = "SYMBOL",
                                                   multiVals = "first"))
  
  # Returns NULL if no Entrez identifier is available.
  if (is.na(entrez)) return(NULL)
  
  # Query KEGG for drugs associated with the corresponding human gene.
  url <- paste0("https://rest.kegg.jp/link/drug/hsa:", entrez)
  resp <- httr::GET(url)
  
  # Skips genes for which the KEGG request is unsuccessful.
  if (httr::http_error(resp)) return(NULL)
  txt <- httr::content(resp, "text")
  
  # Returns NULL when the response does not contain tab-delimited records.
  if (!grepl("\t", txt)) return(NULL)
  
  # Parses the KEGG drug–gene associations into a data frame.
  df <- tryCatch(read.table(text = txt, sep = "\t", stringsAsFactors = FALSE),
                 error = function(e) NULL)
  
  # Skips empty or invalid query results
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  # Standardize the retrieved drug–gene associations and retain corresponding
  # KEGG drug identifiers and source URLs.
  tibble(gene = gene,
         entrez = entrez,
         source_db = "KEGG",
         drug = df$V2,
         source_url = paste0("https://www.kegg.jp/dbget-bin/www_bget?drug:", df$V2))
}


# 23. DrugBank 
#-------------------------------------------------------------------------------
# Extracts drug–target interactions from a locally stored DrugBank XML database 
# and retrieves records corresponding to the requested gene. A valid DrugBank XML 
# database must be available at the path specified inside the function. 
# Processing may be computationally intensive
query_drugbank <- function(gene) {
  
  # The paths for the required data files retrieved from DrugBank stored locally
  xml_file <- "D:/DRYLAB/ATOPICDERMATITIS/Durg Gene-Interaction/full database.xml"
  
  # Parses the DrugBank XML database and extracts drug–target relationships.
  parse_drugbank_gene_interactions <- function(xml_file) {
    
    message("Loading DrugBank XML...")
    doc <- read_xml(xml_file)
    
    # Identifies the XML namespace required for querying DrugBank entries
    ns <- xml_ns(doc)
    
    message("Finding drug entries...")
    drugs <- xml_find_all(doc, "//d1:drug", ns = ns)
    
    message("Extracting targets for ", length(drugs), " drugs...")
    
    # Extracts target information for each DrugBank drug entry.
    interactions <- map_dfr(drugs, function(drug) {
      
      # Retrieves the primary DrugBank identifier and drug name.
      drug_id   <- xml_text(xml_find_first(drug, ".//d1:drugbank-id[@primary='true']", ns = ns))
      drug_name <- xml_text(xml_find_first(drug, ".//d1:name", ns = ns))
      
      # Identifies molecular targets associated with the drug.
      targets <- xml_find_all(drug, ".//d1:targets/d1:target", ns = ns)
      if (length(targets) == 0) return(NULL)
      
      # Extracts gene and target-level information for each drug target.
      map_dfr(targets, function(tg) {
        
        gene <- xml_text(xml_find_first(tg, ".//d1:gene-name", ns = ns))
        
        # Exclude target records without an annotated gene symbol.
        if (is.na(gene) || gene == "") return(NULL)
        
        # Extract and combine reported target actions.
        action_nodes <- xml_find_all(tg, ".//d1:actions/d1:action", ns = ns)
        action <- if (length(action_nodes) == 0) NA_character_
        else paste(xml_text(action_nodes), collapse = ";")
        
        # Standardize the extracted drug–gene interaction information.
        tibble(drugbank_id   = drug_id,
               drug_name     = drug_name,
               gene_symbol   = gene,
               target_action = action,
               organism      = xml_text(xml_find_first(tg, ".//d1:organism", ns = ns)),
               known_action  = xml_text(xml_find_first(tg, ".//d1:known-action", ns = ns)),
               source_url    = paste0("https://go.drugbank.com/drugs/", drug_id))
      })
    })
    
    message("Done! Extracted ", nrow(interactions), " drug–gene interactions.")
    return(interactions)
  }
  
  # Parses the complete DrugBank interaction dataset and retain records
  # corresponding to the gene of interest.
  db_interactions <- parse_drugbank_gene_interactions(xml_file)
  db_interactions %>% filter(gene_symbol == gene)
}


# 24. HCDT
#-------------------------------------------------------------------------------
# Retrieves drug–gene associations from a locally stored HCDT DRUG-GENE.tsv file.
# The file must contain GENE_SYMBOL, DRUG_NAME, PUBCHEM_CID, Datasource, DRUGURL.
# Returns a data.table
query_hcdt <- function(gene_vector) {
  
  # The paths for the required data files retrieved from HCDT and stored locally
  path <- "D:\\DRYLAB\\ATOPICDERMATITIS\\Durg Gene-Interaction\\DRUG-GENE.tsv"
  
  # Loads .tsv
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  
  # Standardize input gene symbols to uppercase for case-insensitive matching
  gene_vector <- toupper(gene_vector)
  dt[, GENE_SYMBOL := toupper(GENE_SYMBOL)]
  
  # Retains drug–gene associations, drug identifiers, source information, and 
  # corresponding URLs to the genes of interest.
  filtered <- dt[GENE_SYMBOL %in% gene_vector]
  
  # Retain drug identifiers, source information, and corresponding URLs.
  result <- filtered[, .(GENE_SYMBOL,
                         DRUG_NAME,
                         PUBCHEM_CID,
                         Datasource,
                         DRUGURL)]
  
  return(result)
}


# 25. TTD
#-------------------------------------------------------------------------------
# Retrieves drug–target associations from the Therapeutic Target Database (TTD) 
# using locally downloaded target, drug cross-matching, and drug–target mapping 
# files. Three TTD files are required: P1-01-TTD_target_download.txt, 
# P1-03-TTD_crossmatching.txt, P1-07-Drug-TargetMapping.xlsx. Returns a dataframe.
query_TTD<- function(genes) {
  
  # Convert gene input to uppercase
  genes <- toupper(genes)
  
  # Defines TTD input files
  target_file <- "D:\\DRYLAB\\ATOPICDERMATITIS\\Durg Gene-Interaction\\P1-01-TTD_target_download.txt"
  drug_file   <- "D:\\DRYLAB\\ATOPICDERMATITIS\\Durg Gene-Interaction\\P1-03-TTD_crossmatching.txt"
  matrix_file <- "D:\\DRYLAB\\ATOPICDERMATITIS\\Durg Gene-Interaction\\P1-07-Drug-TargetMapping.xlsx"
  
  # Maps genes to TTD target identifiers
  txt1 <- readLines(target_file, warn = FALSE)
  txt1 <- txt1[grepl("\t", txt1)]  
  parts1 <- strsplit(txt1, "\t")
  target_df <- data.table(TARGETID = sapply(parts1, `[`, 1),
                          FIELD    = sapply(parts1, `[`, 2),
                          VALUE    = sapply(parts1, function(x) paste(x[3:length(x)], collapse = "\t")))
  
  # Retains TTD records containing gene-name annotations
  target_map <- target_df[FIELD == "GENENAME", .(TARGETID,
                                                 GENE_SYMBOL = toupper(VALUE))]
  
  # Restricts TTD targets to the genes of interest
  targets_of_interest <- target_map[GENE_SYMBOL %in% genes]
  if (nrow(targets_of_interest) == 0) {
    message("No matching TTD targets found for input genes.")
    return(NULL)
  }
  
  # Parses the TTD drug cross-matching file and retains key drug identifiers, 
  # names, and PubChem compound identifiers.
  txt2 <- readLines(drug_file, warn = FALSE)
  txt2 <- txt2[grepl("\t", txt2)]
  parts2 <- strsplit(txt2, "\t")
  drug_df <- data.table(DRUGID = sapply(parts2, `[`, 1),
                        FIELD  = sapply(parts2, `[`, 2),
                        VALUE  = sapply(parts2, function(x) paste(x[3:length(x)], collapse = "\t")))
  
  # Retains the key drug annotation fields required for downstream integration.
  drug_map <- drug_df[FIELD %in% c("TTDDRUID", "DRUGNAME", "PUBCHCID")]
  drug_map <- dcast(drug_map,
                    DRUGID ~ FIELD,
                    value.var = "VALUE")
  
  # Standardizes drug annotation column names.
  setnames(drug_map,
           c("DRUGID", "TTDDRUID", "DRUGNAME", "PUBCHCID"),
           c("DrugID", "TTDDRUID", "DrugName", "PubChemCID"))
  
  # Load the TTD drug–target mapping matrix linking TTD target identifiers 
  # to corresponding drug identifiers.
  TTDmatrix <- read_xlsx(matrix_file)
  TTDmatrix <- data.table(TTDmatrix)
  
  # Integrate target and drug information by linking genes of interest to 
  # their corresponding drug–target associations
  merged1 <- merge(targets_of_interest, TTDmatrix,
                   by.x = "TARGETID", by.y = "TargetID",
                   all.x = TRUE)
  
  # Adds drug names and compound identifiers to the drug–target associations.
  merged2 <- merge(merged1, drug_map,
                   by = "DrugID",
                   all.x = TRUE)
  
  # Generate final TTD interaction dataset. Retains gene, target, drug, 
  # clinical-status, and mechanism-of-action information for downstream 
  # drug–gene interaction analysis.
  result <- merged2[, .(Gene               = GENE_SYMBOL,
                        TargetID           = TARGETID,
                        DrugID             = DrugID,
                        DrugName           = DrugName,
                        PubChemCID         = PubChemCID,
                        Highest_status     = Highest_status,
                        MOA                = MOA)]
  
  return(result)
}



# 26. Drug–gene interaction extraction pipeline 
#-------------------------------------------------------------------------------
# Runs the complete multi-database drug–gene interaction extraction pipeline for 
# a supplied set of genes. It standardizes gene symbols, queries all configured 
# databases, and organizes the results into a nested database-by-gene structure.
dgInt <- function(genes) {
  
  # Standardizes gene symbols and remove duplicate inputs.
  genes <- unique(toupper(genes))
  
  # Runs a single-gene query
  safe_query <- function(fn, gene) {
    tryCatch(
      fn(gene),
      error = function(e) {
        message("[ERROR] ", gene, " in ", deparse(substitute(fn)), ": ", conditionMessage(e))
        tibble()   
      }
    )
  }
  
  # Applies a database-specific query function to each gene and store the
  # resulting interaction records
  build_db_list <- function(db_name, fn) {
    results <- setNames(
      lapply(genes, function(g) safe_query(fn, g)),
      genes
    )
    return(results)
  }
  
  # Builds the nested list
  nested_output <- list(
    DrugCentral  = build_db_list("DrugCentral",query_drugcentral),
    ChEMBL       = build_db_list("ChEMBL",query_chembl),
    PharmGKB     = build_db_list("PharmGKB", query_pharmGKB),
    OpenTargets  = build_db_list("OpenTargets",query_opentargets),
    DrugBank     = build_db_list("DrugBank",query_drugbank),
    KEGG         = build_db_list("KEGG",query_kegg),
    HCDT         = build_db_list("HCDT", query_hcdt),
    TTD         = build_db_list("TTD", query_TTD))
  
  return(nested_output)
}


# 27. Harmonization of drug–gene interaction results 
#-------------------------------------------------------------------------------
# Combines heterogeneous drug–gene interaction results from multiple databases 
# into a standardized master dataset suitable for downstream analysis and 
# prioritization.
# Drugs_out =	Nested multi-database interaction results
# DGIdbres =	Optional DGIdb result table
# save_file =	Optional .xlsx output path
# verbose	= Logical controlling verbose behaviour
# Returns a list containing: Drug_dfs = individula database outputs; 
# Drugs = High scoring interactions from DGidb, DrugCentral and Chembl
# NADrugs = interactions having NA
# Master_Drugs = unified dataframe with results from all databases
harmonize_drug_results <- function( Drugs_out,
                                    DGIdbres = NULL,
                                    save_file = NULL,
                                    verbose = TRUE) {
  # Loads required packages
  library(dplyr)
  library(purrr)
  
  
  # Combine gene-level results from each database while retaining the 
  # originating database as the interaction source.
  Drug_dfs <- imap(Drugs_out, function(db_list, db_name) {
    db_list <- purrr::compact(db_list)   # remove NULLs
    db_list <- db_list[vapply(db_list,function(x) is.data.frame(x) && ncol(x) > 0,logical(1))]
    df <- bind_rows(db_list)
    df$source <- db_name
    # Uses the original HCDT data source when available.
    if (db_name == "HCDT" && "Datasource" %in% colnames(df)) {
      df$source <- df$Datasource
    }
    df
  })
  
  # Map database-specific column names to common fields to enable integration of 
  # heterogeneous drug–gene interaction records.
  Drugs <- lapply(Drug_dfs, function(x) {
    g      <- c("gene_name","gene","gene_symbol","GENE_SYMBOL","Gene")
    drg    <- c("drug_name","drug","DRUG_NAME","DrugName","Drug")
    src    <- c("source","source_db","Datasource","Source")
    int_sc <- c("interaction_score","notes","Interaction score")
    app    <- c("approved","Regulatory approval","Highest_status")
    data.frame(Gene  = as.character(pick_col(x, g)),
               Drug  = as.character(pick_col(x, drg)),
               Source = as.character(pick_col(x, src)),
               Interaction_Score = suppressWarnings(as.numeric(gsub("value: ", "", pick_col(x, int_sc)))),
               Regulatory_Approval_Status = as.factor(pick_col(x, app)), stringsAsFactors = FALSE)      })
  Drugs <- bind_rows(Drugs)
  
  # Separate records without an available interaction score.
  NADrugs <- Drugs[is.na(Drugs$Interaction_Score), ]
  
  # Retain scored interactions and reshape scores by database source.
  Drugs_scored <- Drugs[!is.na(Drugs$Interaction_Score), ]
  Drugs_scored <- Drugs_scored %>%  pivot_wider(., names_from =Source, values_from =Interaction_Score)
  str(Drugs_scored)
  
  # Combines the original database-specific records while retaining the 
  # database identifier for traceability.
  
  Drugs_dfs_merged <- bind_rows(Drug_dfs, .id = "Database")
  
  # Removes database-specific identifiers and metadata that are not required 
  # in the harmonized master table.
  cols_to_remove <- c(
    "interaction_source_db_version",
    "gene_concept_id",
    "drug_concept_id",
    "evidence_sources",
    "pubmed_ids",
    "source_url",
    "evidenceLevel",
    "phenotypeCategory",
    "pmid_count",
    "annotation_id",
    "url",
    "ensembl",
    "drug_id",
    "disease_id",
    "urls",
    "drugbank_id",
    "organism",
    "known_action",
    "entrez",
    "PUBCHEM_CID",
    "DRUGURL",
    "TargetID",
    "DrugID",
    "PubChemCID"
  )
  
  Drugs_dfs_merged <- Drugs_dfs_merged %>%
    dplyr::select(-any_of(cols_to_remove))
  
  # Retrieve a column and return missing values when the column 
  # is not present in a particular database-derived dataset.
  safe_col <- function(df, col) {
    if (col %in% names(df)) {
      df[[col]]
    } else {
      rep(NA_character_, nrow(df))
    }
  }
  
  # Harmonize gene, drug, source, regulatory status, interaction mechanism, 
  # and interaction-score fields across databases to construct the master drug–gene table 
  print("Column names of Drugs_dfs_merged:\n")
  print(str(Drugs_dfs_merged))
  tmp <- Drugs_dfs_merged %>% 
    mutate(across(everything(), as.character)) %>%
    mutate(
      
      # Standardizes database/source information.
      Source = coalesce(
        safe_col(cur_data(), "Datasource"),
        safe_col(cur_data(), "Database"),
        safe_col(cur_data(), "source_db"),
        safe_col(cur_data(), "interaction_source_db_name")
      ),
      
      # Standardizes gene identifiers/symbols.
      Gene = coalesce(
        safe_col(cur_data(), "Gene"),
        safe_col(cur_data(), "GENE_SYMBOL"),
        safe_col(cur_data(), "gene_symbol"),
        safe_col(cur_data(), "gene"),
        safe_col(cur_data(), "gene_claim_name"),
        safe_col(cur_data(), "gene_name")
      ),
      
      # Standardizes drug names.
      Drug = coalesce(
        safe_col(cur_data(), "DrugName"),
        safe_col(cur_data(), "DRUG_NAME"),
        safe_col(cur_data(), "drug_name"),
        safe_col(cur_data(), "drug"),
        safe_col(cur_data(), "drug_claim_name")
      ),
      
      # Standardizes regulatory or clinical development status.
      Regulatory_status = coalesce(
        safe_col(cur_data(), "Highest_status"),
        safe_col(cur_data(), "approved"),
        safe_col(cur_data(), "status"),
        safe_col(cur_data(), "phase"),
        safe_col(cur_data(), "clinical_stage")
      ),
      
      # Standardizes drug–target interaction or mechanism annotations.
      Interaction = coalesce(
        safe_col(cur_data(),"mechanism"),
        safe_col(cur_data(),"target_action"),
        safe_col(cur_data(),"MOA"),
        safe_col(cur_data(),"interaction_type"),
        safe_col(cur_data(),"interaction_types")
      ),
      
      # Retains chemical and evidence-based interaction scores.
      Chemical_Interaction_Score = notes,
      Evidence_Interaction_Score = interaction_score
    )
  
  names(tmp)
  
  # Adds DGIdb interaction scores by matching records using gene and drug names.
  if (!is.null(DGIdbres)) {
    
    Master_Drugs <- tmp %>%
      left_join(
        DGIdbres %>%
          dplyr::rename(
            Gene = gene,
            Drug = drug
          ) %>%
          dplyr::select(
            Gene,
            Drug,
            DGIdb_Interaction_Score
          ),
        by = c("Gene", "Drug")
      )
  }
  
  # Converts interaction scores to numeric values for downstream analysis.
  Master_Drugs <- Master_Drugs %>%
    mutate(
      Evidence_Interaction_Score = suppressWarnings(as.numeric(Evidence_Interaction_Score)),
      
      DGIdb_Interaction_Score = suppressWarnings(as.numeric(DGIdb_Interaction_Score)),
      
      Chemical_Interaction_Score = suppressWarnings(as.numeric(gsub("value:", "", Chemical_Interaction_Score))))
  
  # Replaces missing interaction and regulatory-status annotations with 
  # "unknown" and standardize binary regulatory-status labels.
  Master_Drugs <- Master_Drugs %>%
    mutate( across( c(Interaction, Regulatory_status),~coalesce(., "unknown")),
            Regulatory_status = recode( Regulatory_status, "TRUE"  = "Approved", "FALSE" = "Unapproved"))
  
  
  # Organizes records by gene–drug pair and retain the standardized fields 
  # required for downstream drug prioritization and interpretation.
  Master_Drugs <- Master_Drugs %>%
    group_by(Gene, Drug)
  print(colnames(Master_Drugs))
  str(Master_Drugs)
  Master_Drugs <- Master_Drugs %>%
    dplyr::select(
      Gene,
      Drug,
      Source,
      Regulatory_status,
      Interaction,
      Chemical_Interaction_Score,
      Evidence_Interaction_Score,
      DGIdb_Interaction_Score
    )
  
  # Exports the harmonized master table when an output path is provided.
  if (!is.null(save_file)) {
    writexl::write_xlsx(Master_Drugs, save_file)
  }
  
  # Returns the intermediate and final harmonized datasets.
  return(list(
    Drug_dfs = Drug_dfs,
    Drugs = Drugs_scored,
    NADrugs = NADrugs,
    Master_Drugs = Master_Drugs
  ))
}


# 28.Drug evidence summary
#-------------------------------------------------------------------------------
# Separates approved drug–gene interactions from the complete interaction dataset
# and calculates a standardized chemical binding score from available interaction 
# values. 
# The input data frame must contain: Regulatory_status, Chemical_Interaction_Score
# Chemical interaction values are transformed using:-log10(Chemical_Interaction_Score * 1e-9)
# Returns and saves two excel sheets: Approved_Drug_Gene_Interaction.xlsx, 
# Drug_Gene_Binding_Score.xlsx
drug_evidence_summary2 <- function(Drugs, save_prefix=NULL){
  
  # Retains interactions annotated as approved or with an equivalent TRUE 
  # regulatory-status designation.
  Drugs_app <- Drugs %>%
    dplyr::filter(Regulatory_status %in% c("Approved","TRUE"))
  
  # Retains interactions with available chemical interaction values and 
  # transforms the reported values to a -log10 scale for score normalization.
  Binding_score <- Drugs %>%
    dplyr::filter(!is.na(Chemical_Interaction_Score)) %>%
    dplyr::mutate(Log10ChEMBLintscore = -log10(Chemical_Interaction_Score * 1e-9))
  list( ApprovedDrugs = Drugs_app,
        BindingScore = Binding_score)
  
  # Export approved drug–gene interactions and calculated binding scores for  
  # downstream analysis.
  write_xlsx(Drugs_app, paste0(save_prefix, "_Approved_Drug_Gene_Interaction.xlsx"))
  write_xlsx(Binding_score, paste0(save_prefix, "_Drug_Gene_Binding_Score.xlsx"))
  
  return(list(ApprovedDrugs = Drugs_app,
              BindingScore = Binding_score))
  
}



# 29.Binding score distribution and threshold analysis 
#-------------------------------------------------------------------------------
# Characterizes the distribution of chemical interaction/binding scores, 
# calculates percentile-based thresholds, determines the number of interactions 
# above each threshold, and optionally generates a distribution plot.
# The input data frame must contain a numeric score column Log10ChEMBLintscore
# Binding_score =	Data frame containing interaction scores
# score_col	= Name of score column
# probs =	Percentiles used to calculate thresholds
# make_plot =	Logical; whether to generate ggplot visualization
analyze_binding_scores <- function(Binding_score,
                                   score_col = "Log10ChEMBLintscore",
                                   probs = c(0.25, 0.50, 0.75, 0.95, 0.97, 0.99),
                                   make_plot = TRUE) {
  
  # Loads necessary packages
  library(ggplot2)
  
  # Filters finite and positive interaction scores for downstream distribution
  # and threshold analysis.
  Drugs <- Binding_score[is.finite(Binding_score[[score_col]]) & Binding_score[[score_col]] > 0, ]
  intscore <- Drugs[[score_col]]
  
  # Estimates the score distribution by empirical density of the interaction scores.
  den <- density(na.omit(intscore))
  
  # Calculates predefined score percentiles and determines the number of 
  # interactions meeting or exceeding each corresponding threshold.
  q <- quantile(na.omit(intscore),probs = probs )
  n <- sapply(q, function(thresh) sum(intscore >= thresh, na.rm = TRUE))
  cat("\n=====================================\n")
  cat("Interaction Score Threshold Summary\n")
  cat("=====================================\n")
  for(i in seq_along(q)){
    cat( "\n****** Threshold :",names(q)[i],"percentile\n")
    cat( "****** intscore  :", round(q[i], 6),"\n")
    cat( "****** n         :",n[i],"\n")
  }
  
  # Assess the log-transformed score distribution. Applies a log10 transformation 
  # to the interaction scores and calculates corresponding percentile 
  # thresholds and interaction counts.
  log_scores <- log10(intscore)
  qlog <- quantile(na.omit(log_scores),probs = probs )
  nlog <- sapply(qlog, function(thresh) sum(log_scores >= thresh, na.rm = TRUE))
  cat("\n=====================================\n")
  cat("Log10 Interaction Score Summary\n")
  cat("=====================================\n")
  for(i in seq_along(qlog)) {
    cat("\n****** Threshold :",names(qlog)[i],"percentile\n")
    cat("****** log10(score) :",round(qlog[i], 6),"\n")
    cat("****** n            :",nlog[i],"\n")
  }
  
  # Generates interaction-score density plot
  p <- list()
  plot( den, lwd = 2, main = "Density Plot of Interaction Scores", xlab = "Interaction Scores")
  p$density <- den
  if(make_plot){
    p$Distribution <- ggplot(Drugs, aes(x = Log10ChEMBLintscore)) +
      
      geom_histogram(aes(y = after_stat(density)),
                     bins = 50,
                     fill = "black",
                     color = "gray44",
                     linewidth = 0.6) +
      
      geom_density(color = "#1AAFBC", linewidth = 2, adjust = 0.3) +
      
      geom_vline(xintercept = q[1], linetype = "dashed", linewidth = 1, col = "red") +
      geom_vline(xintercept = q[2], linetype = "dashed", linewidth = 1, col = "red") +
      geom_vline(xintercept = q[3], linetype = "dashed", linewidth = 1, col = "red") +
      geom_vline(xintercept = q[4], linetype = "dashed", linewidth = 1, col = "red") +
      
      annotate("text", x = q[1], y = Inf,
               label = paste0("25% = ", round(q[1]), "\n(n=", n[[1]], ")"),
               vjust = 2, hjust = -0.1, size = 4, col = "turquoise3", fontface = "bold") +
      annotate("text", x = q[2], y = Inf,
               label = paste0("50% = ", round(q[2]), "\n(n=", n[[2]], ")"),
               vjust = 4, hjust = -0.1, size = 4, col = "turquoise3", fontface = "bold") +
      annotate("text", x =q[3], y = Inf,
               label = paste0("75% = ", round(q[3]), "\n(n=", n[[3]], ")"),
               vjust = 6, hjust = -0.1, size = 4, col = "turquoise3", fontface = "bold") +
      annotate("text", x = q[4], y = Inf,
               label = paste0("95% = ", round(q[4]), "\n(n=", n[[4]], ")"),
               vjust = 4, hjust = -0.1, size = 4, col = "turquoise3", fontface = "bold") +
      
      labs(title = "Distribution of Interaction Scores",
           x = "Log10 Interaction Score",
           y = "Density") +
      
      theme_classic() +
      
      theme(plot.title = element_text(hjust =0.5, face = "bold", size = 16),
            axis.title.x = element_text(hjust =0.5, face = "bold", size = 12),
            axis.title.y = element_text(hjust =0.5, face = "bold", size = 12)) +
      print(p$Distribution)
  }
  
  # Returns filtered scores, density estimates, percentile thresholds, 
  # interaction counts, and generated plots.
  return(list(Cleaned_Drugs = Drugs,
              Density = den,
              Quantiles = q,
              AbsCounts = n,
              LogQuantiles = qlog,
              LogCounts = nlog,
              Plot = p))
}




