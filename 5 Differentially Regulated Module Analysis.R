#===============================================================================
#                *** Differentially Regulated Module Analysis ***            
#-------------------------------------------------------------------------------
#                                *** WGCNA ***                     
#===============================================================================


# **** Packages
{
  library(WGCNA)
  library(ggplot2)
  library(tibble)
  library(limma)
  library(dplyr)
  library(gridExtra)
  library(tidyverse)
  library(tidyr)
  library(ggforce)
  library(rlang)
  library(broom)
  library(igraph)
  library(ggraph)
  library(tidyr)
  library(openxlsx)
  library(FSA)
  library(effsize)
  library(writexl)
  library(stringr)
}



# **** WGCNA ****
#=============================================================================== 



# 1. Data Preprocessing and quality control
#-------------------------------------------------------------------------------
{
  # Loading the batch corrected master expression dataframe and transpose it to 
  # obtain samples as rows and genes as columnsa as the input orientation required for WGCNA.
  MasterAD <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Network Analysis\\Master_Dataset.rds")
  tMasterAD <- t(MasterAD)
  # [1]   826 13859
}
  
  
  
# 2. Quality control and sample outlier detection 
#-------------------------------------------------------------------------------
{ 
  # Assessing genes and samples quality using WGCNA's goodSamplesGenes() function.
  gsg <- goodSamplesGenes(tMasterAD)
  summary(gsg)
  #             Length Class  Mode   
  # goodGenes   13859  -none- logical
  # goodSamples   826  -none- logical
  # allOK           1  -none- logical
  # gsg$allOK
  
  table(gsg$goodGenes)
  # TRUE 
  # 13859 
  
  table(gsg$goodSamples)
  # TRUE 
  # 826
  
  # PCA-based assessment to visualize sample-level structure and identify potential 
  # transcriptional outliers
  pca <- prcomp(t(filMasterAD))
  pca.dat <- pca$x
  pca.var <- pca$sdev^2
  pca.var.percent <- round(pca.var/sum(pca.var)*100, 
                           digits = 2)
  pca.dat <- as.data.frame(pca.dat)
  ggplot(pca.dat, aes(PC1, PC2)) +
    geom_point() +
    geom_text(label = rownames(pca.dat)) +
    labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
         y = paste0('PC2: ', pca.var.percent[2], ' %'))
  
  
  # Removing samples identified as outliers based on the PCA criteria plot shared in Fig. S2a.
  exsamples <- rownames(pca.dat[pca.dat$PC1 < -1.3e6 | pca.dat$PC2 < -1000000 | pca.dat$PC2 > 625000, ])
  # [1] "GSM4874250" "GSM5788596" "GSM5788624" "GSM5788862" "GSM8084719"
  # [6] "GSM8535188" "GSM8535190" "GSM8535191" "GSM8535200"
  filMasterAD <- filMasterAD[,!(colnames(filMasterAD) %in% exsamples)]         
  dim(filMasterAD)
  # [1] 11052   817
  } 
  
  
  
# 3. Expression filtering and variance normalization 
#-------------------------------------------------------------------------------
{
  # Genes are retained if their expression reaches the minimum count threshold 
  # in at least 75% of the samples. Removing lowly expressed genes to reduce noise 
  # and improve the robustness of correlation-based network construction. 
  filMasterAD <- filMasterAD[rowSums(filMasterAD >= 15) >= 612.75,]
  nrow(filMasterAD)
  # 11052 genes

  # z-score transformation of expression profiles across samples, resulting in 
  # mean-centered and unit-variance expression profiles for network analysis.
  # z-score per gene (mean=0, sd=1)
    filMasterAD <- scale(t(filMasterAD))
    saveRDS(filMasterAD, "FilteredDataset.rds")
    
  # redefining objects
    TfilMasterAD <- filMasterAD
    filMasterAD <- t(filMasterAD)
  }
  
  
  
# 4. Weighted gene co-expression network construction Network Construction  
#-------------------------------------------------------------------------------
{
  # A set of soft-thresholding powers
  power <- c(1:30)
  
  # Evaluating candidate soft-thresholding powers based on their scale-free
  # topology fit and corresponding network connectivity and extracting the fit.
  sft <- pickSoftThreshold(TfilMasterAD,
                           powerVector = power,
                           networkType = "signed",
                           verbose = 5)
  sft.data <- sft$fitIndices
  
  
  # Selecting the soft-thresholding power based on the scale-free topology fit and
  # mean network connectivity.
  # Visualization of scale-free topology fit and mean network connectivity across to pick power
  
  # R^2 plot
  a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
        geom_point() +
        geom_text(nudge_y = 0.1) +
        geom_hline(yintercept = 0.8, color = 'red') +
        labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
        theme_classic()
  
  # k-mean plot
  a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
        geom_point() +
        geom_text(nudge_y = 0.1) +
        labs(x = 'Power', y = 'Mean Connectivity') +
        theme_classic()
  
  # Arrange the scale-free topology fit and mean connectivity plots vertically Fig.S2b
  grid.arrange(a1, a2, nrow = 2)
  
  # Chosen power
  soft_power <- 20
  
  
  # Constructing a signed weighted gene co-expression network using the selected
  # soft-thresholding power.Fix the random seed and disable parallel computation 
  # to improve reproducibility of network construction.
  set.seed(1234)                                             # Sets RNG seed
  disableWGCNAThreads()                                      # Disables multithreading (forces deterministic order)
  temp_cor <- cor                                            # Temporarily assigns WGCNA's correlation function to ensure compatibility
  cor <- WGCNA::cor                                          # with downstream network construction and topological overlap calculations.                 
  bwnet <- blockwiseModules( TfilMasterAD,                   # Running blockwiseModules with fixed parameters
                             power = soft_power,
                             TOMType = "signed",
                             networkType = "signed",
                             maxBlockSize = 14000,
                             deepSplit = 2,
                             minModuleSize = 30,
                             mergeCutHeight = 0.25,
                             numericLabels = FALSE,
                             pamRespectsDendro = TRUE,
                             saveTOMs = TRUE,
                             saveTOMFileBase = "TOM",
                             randomSeed = 1234,
                             verbose = 3)
  cor <- temp_cor                                            # Restores the original correlation function after completing the WGCNA calculations
  View(bwnet)
  
  
  # Calculating the topological overlap matrix (TOM), which incorporates both direct
  # gene-gene correlations and shared network connectivity.
  set.seed(1234)                                             # Sets RNG seed
  disableWGCNAThreads()                                      # Disables multithreading (forces deterministic order) or allowWGCNAThreads(1)
  Sys.setenv(OMP_NUM_THREADS = 1)                            # Disables system-level multithreading for BLAS/LAPACK
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)                       # (prevents OpenBLAS/MKL from parallelizing matrix ops)
  Sys.setenv(MKL_NUM_THREADS = 1)
  Sys.setenv(NUMEXPR_NUM_THREADS = 1)
  TOM <- TOMsimilarityFromExpr( datExpr = TfilMasterAD,      # Generating topological overlap matrix
                                weights = NULL,
                                corType = "pearson",
                                networkType = "signed",
                                power = soft_power,
                                TOMType = "signed",
                                TOMDenom = "min",
                                suppressNegativeTOM = FALSE,
                                useInternalMatrixAlgebra = FALSE,
                                verbose = 5,
                                indent = 1)
  
  colnames(TOM) <- rownames(filMasterAD)
  rownames(TOM) <- rownames(filMasterAD)
  identical(rownames(TOM), colnames(TOM))
  # [1] TRUE
  # Adjacency matirx
  adjacency <- adjacency(TfilMasterAD, power = soft_power, type = "signed")
  
  
  # Gene dendrogram vizualization  with module assignments before and
  # after module merging to assess the resulting module structure.  
  plotDendroAndColors(bwnet$dendrograms[[1]], 
                      cbind(bwnet$unmergedColors, 
                            bwnet$colors),
                      c("unmerged", "merged"),
                      dendroLabels = FALSE,
                      addGuide = TRUE,
                      hang= 0.03,
                      guideHang = 0.05)
  
  }
  
  
  
# 5. Module eigengene analysis
#-------------------------------------------------------------------------------
{
  module_eigengenes <- bwnet$MEs
  head(module_eigengenes)
  table(bwnet$colors)
  #  black         blue        brown         cyan        green 
  #    277          768          644           49          323 
  # greenyellow         grey    lightcyan      magenta midnightblue 
  #          78         5422           32          147           47 
  # pink       purple          red       salmon          tan 
  #  241          145          308           56           66 
  # turquoise       yellow 
  #      1829          620
  modules <- sort(unique(bwnet$colors))
  
  # Module eigengene dendrogram
  # Calculating pairwise module eigengene dissimilarity as one minus the Pearson
  # correlation between module eigengenes.Performing hierarchical clustering of 
  # module eigengenes.
  MEDiss <- 1 - cor(module_eigengenes)
  METree <- hclust(as.dist(MEDiss), method = "average")
  plot( METree,  main = "Clustering of module eigengenes",  xlab = "", sub = "")
  
  
  # Association between module eigengenes and disease condition.
  # Healthy samples are specified as the reference level for the condition variable.
  metadata <- readRDS("D:/DRYLAB/ATOPICDERMATITIS/Network Analysis/masterMetadata.rds")
  metadata <- metadata %>% filter(!geo_accession %in% exsamples)
  metadata$Condition <- relevel(metadata$Condition, ref = "Healthy")
  head(metadata)
  # Dataset geo_accession Condition
  # GSM2741987 GSE102628    GSM2741987  Lesional
  # GSM2741988 GSE102628    GSM2741988  Lesional
  # GSM2741989 GSE102628    GSM2741989  Lesional
  # GSM2741990 GSE102628    GSM2741990  Lesional
  # GSM2741991 GSE102628    GSM2741991  Lesional
  # GSM2741992 GSE102628    GSM2741992  Lesional
  all.equal(metadata$geo_accession, rownames(module_eigengenes))
  # [1] TRUE
  
  
  # Fitting a linear model to test differential module eigengene expression across
  # skin conditions and apply empirical Bayes moderation to improve estimation
  # of standard errors.
  des_mat <- model.matrix(~Condition, data = metadata[!metadata$geo_accession %in% exsamples, ])
  fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
  fit <- limma::eBayes(fit)
  
  
  # Apply multiple-testing correction and extract statistical summaries for
  # module-condition associations.
  stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
              tibble::rownames_to_column("module")
  write.xlsx(stats_df, "Module-trait relationship.xlsx", overwrite = TRUE)
    
  
  # Reshaping module eigengene data into long format for visualization across
  # disease conditions.  
  modules_meta <- module_eigengenes %>% 
                  tibble::rownames_to_column("geo_accession" )%>%
                  inner_join(metadata %>% 
                               select(geo_accession, Condition), 
                             by = "geo_accession")
  modules_long <- modules_meta %>%
                  pivot_longer(cols = starts_with("ME"),       
                               names_to = "Module",           
                               values_to = "Eigengene")
  
  
  # Boxplot of distribution of module eigenegenes across skin conditions 
  # with adjusted.(Fig. 3a)
  ggplot(modules_long, aes(x = Condition, y = Eigengene,fill = Condition)) +
        geom_boxplot(width = 0.5, outlier.shape = NA) +
        facet_wrap(~ Module, scales = "free_y", ncol = 4) +
        geom_text( data = stats_df %>%
                          group_by(module) %>%
                          reframe(pvalue = unique(adj.P.Val), .groups = "drop") %>%
                          dplyr::rename(Module = module), 
                   aes(label = paste0("p = ", signif(pvalue, 3))),
                   x = Inf, y = Inf,
                   hjust = 1.1, vjust = 1.5,
                   inherit.aes = FALSE, size = 8) +
        theme_classic() +
        labs(title ="Distribution of Module Eigengene in Skin Conditions",
             y = "Module Eigengene") +
        theme(
          plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size =16, face = "bold"),
          axis.text.y = element_text(size = 18),
          strip.text = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16, face = "bold"), 
          legend.title = element_text(size = 18, face = "bold")
                                                                )
  }
 


# 6. Modules associated with skin condition
#-------------------------------------------------------------------------------
{
  # 1.) Modules associated with skin conditions : Kruskal- Wallis Test
  {
    # Applying Kruskal-Wallis test to every module
    k.test <- lapply(names(module_eigengenes), function(mod){
      x <- module_eigengenes[[mod]]
      test <- kruskal.test(x ~ metadata$Condition)
      Kruskal_res <- data.frame(Module = mod,
                                Statistic = as.numeric(test$statistic),
                                Pvalue = test$p.value)
      # Benjaminin Hochberg correction of p values
      Kruskal_res$P.adj <- p.adjust(Kruskal_res$Pvalue, method = "BH")
      return(Kruskal_res)
    })
    Kruskal_res <- do.call(rbind, k.test) 
    
    # Removing the "ME" prefix
    Kruskal_res$Module <- gsub("^ME", "", Kruskal_res$Module)
    
    # Retaining modules showing statistically significant differences after 
    # multiple-testing correction.
    Kruskal_res <- Kruskal_res %>% filter(P.adj <0.01)
    write.xlsx(Kruskal_res, "Kruskal-Wallis test.xlsx", overwrite = TRUE)
    
    # Bar plot of Kruskal-Wallis test statistics
    ggplot(Kruskal_res, aes(x = Module, y = Statistic, fill = Module)) +
      geom_col() +
      scale_fill_manual(values = module_colors) +
      labs(title = "Kruskal–Wallis Statistic pvalue < 0.01",
           x = "Module",
           y = "KW Statistic") +
      theme_classic() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 7)) +
      theme( axis.title.x = element_text(face = "bold", size = 18),
             axis.title.y = element_text(face = "bold", size = 18),
             plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
             axis.text.x = element_text(hjust = 0.5, size = 8),
             axis.text.y = element_text(hjust = 1, size = 16),
             legend.position = "none")
    
    write_xlsx(Kruskal_res, "Module Kruskal-Wallis Test.xlsx")
  }
  
  
  
  # 2.) Modules associated with skin conditions : pairwise post-hoc Dunn Test
  {
    # Defining the order of sample conditions for consistent pairwise comparisons.
    metadata$Condition <- factor( metadata$Condition, levels = c("Healthy", "Lesional", "Non-lesional"))
    
    # Applying pairwise Dunn tests following the Kruskal–Wallis analysis to every module
    d.test <- lapply(names(module_eigengenes), function(mod){
      
      x <- module_eigengenes[[mod]]
      
      # Dunn test
      dtest <- dunnTest(x ~ metadata$Condition)
      dres <- dtest$res
      
      # Adding module name to each row
      dres$Module <- mod
      
      # Benjaminin Hochberg correction of p values
      dres$P.adj <- p.adjust(dres$P.unadj, method = "BH")
      
      return(dres)
    })
    Dunn_res <- do.call(rbind, d.test)
    
    # Standardizing comparison pairs and orienting Z statistics relative to the specified biological comparisons.
    Dunn_res <- Dunn_res %>%
      mutate( Comparison = case_when( Comparison == "Healthy - Lesional" ~ "Lesional vs Healthy",
                                      Comparison == "Healthy - Non-lesional" ~ "Non-lesional vs Healthy",
                                      TRUE ~ Comparison ),
              Z = ifelse( Comparison %in% c("Lesional vs Healthy", "Non-lesional vs Healthy"), -Z, Z )) %>%
      filter( Comparison %in% c( "Lesional vs Healthy", "Non-lesional vs Healthy", "Lesional - Non-lesional")) %>%
      mutate( Comparison = recode( Comparison, "Lesional - Non-lesional" = "Lesional vs Non-lesional"))
    
    # Retaining modules showing statistically significant differences after
    # multiple-testing correction.
    Dunn_res <- Dunn_res %>% filter(P.adj <0.05)
    
    # Removing the "ME" prefix
    Dunn_res <- Dunn_res %>% mutate(Module = gsub("ME", "", Module))
    write_xlsx(Dunn_res, "Module Dunn Test.xlsx")
    
    # Bar plot of Kruskal-Wallis test statistics
    ggplot(Dunn_res, aes(x = Module, y = Z, fill = Module)) +
      geom_col( colour = "black", linewidth = 0.3) +
      scale_fill_manual(values = module_colors) +
      labs(  title = "Dunn Test pvalue < 0.05",
             x = "Module",
             y = "Dunn test Z score") +
      theme_classic() +
      theme( axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
             axis.text.y = element_text(size = 12),
             axis.title.x = element_text(face = "bold", size = 14),
             axis.title.y = element_text(face = "bold", size = 14),
             plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
             legend.position = "none",
             strip.text =  element_text(face = "bold", size = 12),
             strip.background = element_rect(fill = "#F2F2F2") ) +
      facet_wrap(~Comparison)
    
    write_xlsx(Dunn_res, "Module Dunn Test.xlsx")
  }
  
  
  
  # 3.) Modules associated with skin conditions : pairwise effect sizes using Cliff's delta
  {
    # Defining the three skin condition pairwise comparisons.
    pairs <- list( c("Lesional", "Healthy"),
                   c("Non-lesional", "Healthy"),
                   c("Lesional", "Non-lesional"))
    
    # Estimating  Cliff's delta for each module across all pairwise comparisons.
    Cliff_res <- lapply(colnames(module_eigengenes), function(mod){
      
      ME  <- module_eigengenes[[mod]]
      grp_all <- metadata$Condition   
      
      res <- lapply(pairs, function(p){
        
        class1 <- p[1]
        class2 <- p[2]
        
        idx <- grp_all %in% p
        x   <- ME[idx]
        g   <- grp_all[idx]
        
        # Removing NA (missing module eigengene or group values)
        keep <- complete.cases(x, g)
        x <- x[keep]
        g <- g[keep]
        
        # Defining the comparison order: class1 > class2
        g <- factor(g, levels = c(class1, class2))
        
        # Returning an undefined effect size if both groups are not represented.
        if (length(unique(g)) < 2) {
          return(data.frame( Module = mod,
                             Comparison = paste(class1, "vs", class2),
                             CliffDelta = NA))
        }
        
        # Assigning zero effect size when both groups contain only a single
        # identical value.
        if (length(unique(x[g == class1])) == 1 &&
            length(unique(x[g == class2])) == 1) {
          return(data.frame(Module = mod,
                            Comparison = paste(class1, "vs", class2),
                            CliffDelta = 0))
        }
        
        # Estimating Cliff's delta .
        cd <- tryCatch(cliff.delta(x ~ g),
                       error = function(e) NULL)
        
      })
      
      do.call(rbind, res)
    })
    Cliff_res <- do.call(rbind, Cliff_res)
    
    # Removing the "ME" prefix
    Cliff_res <- Cliff_res %>% mutate(Module = gsub("ME", "", Module))
    
    # Bar plot of Cliff's Delta
    ggplot(Cliff_res, aes(x = Module, y = CliffDelta, fill = Module)) +
      geom_col( colour = "black", linewidth = 0.3) +
      scale_fill_manual(values = module_colors) +
      labs(  title = "Cliff's Delta",
             x = "Module",
             y = "Cliff's Delta") +
      theme_classic() +
      theme( axis.text.x = element_text(angle = 45, hjust = 1),
             axis.title.x = element_text(face = "bold", size = 14),
             axis.title.y = element_text(face = "bold", size = 14),
             plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
             legend.position = "none",
             strip.text =  element_text(face = "bold", size = 12),
             strip.background = element_rect(fill = "#F2F2F2") ) +
      facet_wrap(~Comparison)
    
    write_xlsx(Cliff_res, "Module Cliffs Delta.xlsx")
    
  }
  
  # ..............................................................................
} 
 


# 7. Intramodular analysis: GS, MM
#-------------------------------------------------------------------------------
{
  # Encoding disease conditions as binary category
  nSamples <- nrow(TfilMasterAD)
  traits <- metadata %>%
            mutate(disease_state_bin = ifelse(Condition == "Lesional", 1,
                                              ifelse(Condition == "Non-lesional", 2, 0)))
  
  # Calculating the gene significance (GS) as the correlation between gene expression and
  # the trait along with associated p-values
  gene.signf.corr <- as.data.frame(cor(TfilMasterAD, traits$disease_state_bin, use = 'p'))
  gene.signf.corr$pvals <- corPvalueStudent(as.matrix(gene.signf.corr), nSamples)[,1]
  colnames(gene.signf.corr) <- c("Coefficient", "pval")
  gene.signf.corr <-  gene.signf.corr %>% rownames_to_column(var = "Gene")
  gene.signf.corr$Module <- bwnet$colors
  
  
  # Calculating module membership (kME or MM) as the correlation between each gene
  # expression and the corresponding module eigengene.
  module.membership.measure <- as.data.frame(t(cor(module_eigengenes, TfilMasterAD, use = 'p')))
  module.membership.measure <- module.membership.measure %>% rownames_to_column(var = "Gene")
  rownames(module.membership.measure) <- module.membership.measure$Gene
  module.membership.measure.pvals <- as.data.frame(corPvalueStudent(as.matrix(module.membership.measure), nSamples))
  
}
  
  
  
# 8. Visualization of intramodular relationships
#-------------------------------------------------------------------------------
{
  
  # Gene Significance (GS)
  #-----------------------
  gs <- gene.signf.corr %>%
    filter(pval<0.05)%>%
    group_by(Module)
  
  # Mean GS
  gene.signf.module.mean <- gs %>%
    summarise(Mean = mean(abs(Coefficient), na.rm = TRUE))
  write.xlsx(gene.signf.module.mean, "Module-trait Gene-Significance.xlsx", overwrite = TRUE)
  par(mar = c(8,10,4,10), mgp = c(6, 1, 0))
  barplot(height = as.matrix(gene.signf.module.mean[, 2]),
          beside = TRUE,
          col = gene.signf.module.mean$Module,
          names.arg = gene.signf.module.mean$Module,
          las = 2,                 # rotate x labels
          ylab = "Mean of Gene Significance Coefficient (p-value < 0.05)",
          ylim = c(0,0.15),
          xlab = "Modules",
          main = "Mean Gene Significance per Module",
          cex.names = 0.8,
          width = 0.1)
  
  # Distribution of GS
  # dynamically color mappings 
  border_colors <- setNames(alpha(modules, 1), modules)   # solid border
  fill_colors   <- setNames(alpha(modules, 0.3), modules) # semi-transparent fill
  par(mar=c(7,4,2,1), mgp = c(3, 0.7, 0)) 
  box <- boxplot(abs(Coefficient) ~ Module, 
                 data = gs,
                 boxwex = 0.6,
                 notch = TRUE,
                 main = "Distribution of Gene Significance per Module ",
                 ylab = "",
                 xlab = "",
                 cex.lab= 1,
                 font.lab = 4,
                 outline = FALSE,
                 las = 2,  
                 col = fill_colors ,   
                 border = border_colors)
  abline(h = max(box$stats[3,]), col = "red")
  stripchart(abs(Coefficient) ~ Module, 
             data = gs,
             method = "jitter",
             pch = 20,
             col = border_colors,   
             vertical = TRUE,
             add = TRUE)
  median_values <- box$stats[3, ]
  mtext("Module", side = 1, line = 5, cex = 1, font = 2)
  mtext( "|Gene Significance| (p-value < 0.05)",
         side = 2,
         line = 3,
         cex = 1,
         font = 2)
  
  gene.signf.module.median <- data.frame(Module = modules,
                                         Median_Gene_Significance = median_values)
  write.xlsx(gene.signf.module.median, "Module-trait Gene-Significance Median.xlsx", overwrite = TRUE)
  
  
  
  
  # GS vs MM relationship
  #-----------------------
    gsmm <- gene.signf.corr %>% select(1:4)
    # Adding MM columns for every gene
    for (mod in names(module_eigengenes)) {
      gsmm[[paste0("MM.", mod)]] <- module.membership.measure[[mod]]
      gsmm[[paste0("pval.", mod)]] <- module.membership.measure.pvals[[mod]]
    }
    gsmm_long <- gsmm %>%
      dplyr::rename(Coefficient.pval = pval) %>%
      select(Gene, Coefficient, Coefficient.pval, Module, matches("^MM\\.ME"), matches("^pval\\.ME")) %>%
      pivot_longer( cols = matches("^(MM|pval)\\.ME"),
                    names_to = c("Type", "ModuleColor"),
                    names_pattern = "(MM|pval)\\.ME(.*)",
                    values_to = "value") %>%
      pivot_wider( names_from = Type,
                   values_from = value) %>%
      dplyr::rename( Module.membership = MM,
                     Module.membership.pval = pval,
                     module = ModuleColor) %>%
      select(Gene, Coefficient, Coefficient.pval, module, Module.membership, Module.membership.pval) %>%
      filter(Module.membership.pval < 0.05 & Coefficient.pval < 0.05)
    
    # Linear regression between gene significance and module membership.
    gsmm_slope <- gsmm_long %>%
      group_by(module) %>%
      do(broom::tidy(lm(Coefficient ~ Module.membership, data = .))) %>%
      filter(term == "Module.membership") %>%
      arrange(desc(estimate))
    write.xlsx(gsmm_slope[,-2], "Module-GS-MM regression.xlsx", overwrite = TRUE)
    
    ggplot(slopes_df,
           aes(x = reorder(module, estimate),
               y = estimate,
               fill = module)) +
      geom_col(width = 0.7, colour = "black") +
      coord_flip() +
      scale_y_continuous(limits = c(-0.10, 0.15),
                         breaks = seq(-0.10, 0.15, by = 0.05)) +
      scale_fill_identity() +
      labs(title = "Relation of Gene-Significance with Module-Membership",
           x = "Modules",
           y = "Slope of GS~MM regression") +
      theme_classic() +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 18, face = "bold"),
            axis.text = element_text(size = 14),
            legend.position = "none")
  
}



# 9. Modules Preservation in skin conditions
#-------------------------------------------------------------------------------
{
# a.) Preparing reference and test expression matrices
{
  options(stringsAsFactors = FALSE)
  
  # The lesional samples as the reference dataset and healthy and
  # non-lesional samples as the test dataset for module preservation analysis.
  refExpr  <- filMasterAD[ , metadata$Condition == "Lesional"]  
  testExpr <- filMasterAD[ , metadata$Condition %in% c("Healthy","Non-lesional")]     
  
  # Retaining genes common to both datasets and ensuring identical gene ordering
  # before comparing network modules across datasets.
  common_genes <- intersect(rownames(refExpr), rownames(testExpr))
  length(common_genes)  # 11052
  refExpr  <- refExpr[common_genes, ]
  testExpr <- testExpr[common_genes, ]
  
  # Transpoinf the expression matrix to obtain samples as rows and genes as
  # columns
  TrefExpr <- t(refExpr)
}



# b.) Evaluating module preservation between reference and test datasets
{
  # Multi-set expression object using the lesional network as the reference and 
  # healthy/non-lesional samples as the test dataset.
  multiExpr <- list( ref  = list(data = TrefExpr),
                     test = list(data = t(testExpr)))
  
  # Assigning modules to the reference network for their evaluation
  colorList <- list(ref = bwnet$colors)
  
  # Check: gene in reference network and WGCNA network.
  all(colnames(TrefExpr) == names(bwnet$colors))
  # [1] TRUE
  
  # Estimating preservation statistics using module preservation() at 2,000 permutations 
  pres2000 <- modulePreservation( multiExpr,
                                  colorList,
                                  referenceNetworks = 1,
                                  nPermutations = 2000,
                                  verbose = 3)
}



# c.) Extraction of preservation statistics:  Z-summary and MedianRank
{
  ref <- 1
  test <- 2
  
  # Extracting Z-summary
  stats2000 <- pres2000$preservation$Z$ref[[test]]
  Zsummary2000 <- stats2000$Zsummary.pres
  
  # Extracting median rank statistics
  rank_stats2000 <- pres2000$preservation$observed$ref[[test]]
  medianRank2000 <- rank_stats2000$medianRank.pres
  
  # Compilation of preservation statistics
  res2000 <- data.frame( Module = rownames(stats2000),
                         Zsummary = Zsummary2000,
                         MedianRank = medianRank2000 )
  
  res2000[order(res2000$Zsummary, decreasing = TRUE), ]
  #          Module  Zsummary MedianRank
  # 17    turquoise 42.089030         14
  # 3         brown 37.510423          6
  # 18       yellow 35.480815          6
  # 2          blue 34.034564         10
  # 1         black 29.097769          1
  # 14          red 27.765894          2
  # 6         green 21.858443          7
  # 12         pink 21.093139         11
  # 10      magenta 19.994823          4
  # 5          gold 17.848738         17
  # 8          grey 17.449703         17
  # 7   greenyellow 12.943100          5
  # 13       purple 12.404089         16
  # 16          tan 11.592301          8
  # 15       salmon  9.200960          9
  # 4          cyan  8.032149         10
  # 11 midnightblue  6.912213         13
  # 9     lightcyan  5.272531         14
  
  writexl::write_xlsx(res2000, "Module Preservation.xlsx")
  
}



# d.) Interpret
{
  subset(res2000, Zsummary < 2)       # Non-preserved modules Strongly differential modules
  # [1] Module     Zsummary   MedianRank
  # <0 rows> (or 0-length row.names)
  
  subset(res2000, Zsummary < 10)      # Weakly preserved modules, Weakly differential modules
  # Module Zsummary MedianRank
  # 4          cyan 8.032149         10
  # 9     lightcyan 5.272531         14
  # 11 midnightblue 6.912213         13
  # 15       salmon 9.200960          9
  
}



# e.) Visualization of module preservation
{
  # Defining consistent colours for visualization of individual WGCNA modules.
  module_colors <- c(  "turquoise"    = "turquoise",
                       "brown"        = "brown",
                       "yellow"       = "yellow",
                       "blue"         = "blue",
                       "black"        = "black",
                       "red"          = "red",
                       "pink"         = "pink",
                       "green"        = "green",
                       "magenta"      = "magenta",
                       "grey"         = "grey",
                       "greenyellow"  = "greenyellow",
                       "purple"       = "purple",
                       "tan"          = "tan",
                       "salmon"       = "salmon",
                       "cyan"         = "cyan",
                       "midnightblue" = "midnightblue",
                       "lightcyan"    = "lightcyan",
                       "gold"         = "gold")
  
  # Scatter plot of module preservation using median rank and Z-summary statistics.
  # Dashed lines indicate the conventional thresholds for weak and strong
  # module preservation.
  ggplot(res2000, aes(x = MedianRank, y = Zsummary, fill = Module)) +
    geom_point(shape = 21, size = 6, color = "black", stroke = 1.2) +
    geom_hline(yintercept = 2, color = "grey", linetype = "dashed", size = 1) +
    geom_hline(yintercept = 10, color = "grey", linetype = "dashed", size = 1) +
    scale_fill_manual(values = module_colors) +
    scale_x_continuous(limits = c(1, 20), breaks = seq(1, 20, by = 1)) +
    scale_y_continuous(breaks = seq(0, max(res200$Zsummary) + 5, by = 5)) +
    labs(title = "Module Preservation",
         x = "Median Rank",
         y = "Z-summary") +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, size = 15, face = "bold"),
          axis.text.y = element_text(hjust = 1, size = 15, face = "bold"),
          axis.title.x = element_text(face = "bold", size = 18),
          axis.title.y = element_text(face = "bold", size = 18),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          legend.position = "right",
          legend.text = element_text(size = 14),   
          legend.title = element_text(size = 18, face = "bold"))
  
}
  
#...............................................................................
}



# 10. Genes in the identified modules
#-------------------------------------------------------------------------------
{
  identified_modules <- c("lightcyan", "midnightblue")
  GeneSymbols <- readRDS("D:/DRYLAB/ATOPICDERMATITIS/Mapping/GeneAnnotation_ensembl.v86_corrected.rds")
  
  module_genes <- lapply(identified_modules, function(x) {
    
    # Extracting and annotating genes of the selected modules
    df <- tibble::enframe(bwnet$colors, name = "Genes", value = "module") %>%
      filter(module == x) %>%
      full_join(GeneSymbols, by = c("Genes" = "ENSEMBL")) %>%
      slice_head(n = rows)
    
    # Assigning LOC identifiers to unannotated genes and removing incomplete entries.
    df$SYMBOL[is.na(df$SYMBOL)] <- paste0("LOC", df$Genes[is.na(df$SYMBOL)])
    df <- na.omit(df)
    rownames(df) <- df$SYMBOL
    
    as.data.frame(df)
    
  })
  
  names(module_genes) <- identified_modules
  
  # Exporting module gene lists to Excel.
  library(openxlsx)
  wb <- createWorkbook()
  for (x in names(module_genes)) {
    addWorksheet(wb, x)
    writeDataTable(wb, x, module_genes[[x]], startRow = 3, startCol = 2)
  }
  saveWorkbook(wb, "Identified Modules.xlsx", overwrite = TRUE)
}


  
# 11. Edge and Node file for selected modules
#-------------------------------------------------------------------------------
{
  # Function to export edge files for gene sets using blockwiseModules_network/
  # TOM matrix 
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
  
  # Function to crossmap the identifiers from Ensembl to Symbols and vice-versa
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
  GeneAnnot <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Mapping\\GeneAnnotation_ensembl.v86_corrected.rds")
  
  # Exporting edge files for  module lightcyan
    cyt.lightcyan <- GenesetCytoscapeExport(geneset = "lightcyan",
                                              TOM = TOM,
                                              blockwiseModules_network = bwnet,
                                              GeneAnnotation = GeneAnnot,
                                              geneset_type = "module",
                                              thresholds = seq(0, 0.1, 0.001),
                                              export = TRUE,
                                              export_threshold = 0.02,
                                              prefix = "Lightcyan",
                                              file_format = ".txt")
    
   # Exporting edge files for module midnightblue
    cyt.midnightblue <- GenesetCytoscapeExport(geneset = "midnightblue",
                                              TOM = TOM,
                                              blockwiseModules_network = bwnet,
                                              GeneAnnotation = GeneAnnot,
                                              geneset_type = "module",
                                              thresholds = seq(0, 0.02 , by = 0.001),
                                              export = TRUE,
                                              export_threshold = 0.016,
                                              prefix = "Midnightblue",
                                              file_format = ".txt")
    
    # Exporting edge files for Lesional vs Healthy DEGs
    metaL <- readRDS("D:/DRYLAB/ATOPICDERMATITIS/Validation/meta-analysis_L_DEG.rds")
    {
      # Upregulated Genes
      inmoduleUPmetaL <- metaL %>% filter(beta >= 1)
      genesDEGmetaLup <- crossmap(inmoduleUPmetaL, id_col = "GeneSymbol", from = "symbol", to = "ensembl")
       # Mapped: 432 / 432 
      genesDEGmetaLup <- intersect(rownames(genesDEGmetaLup), rownames(TOM))
      length(genesDEGmetaLup)
       # [1] 209
      cyt.metaLup <- GenesetCytoscapeExport(geneset = genesDEGmetaLup,
                                            TOM = TOM,
                                            blockwiseModules_network = bwnet,
                                            GeneAnnotation = GeneAnnot,
                                            geneset_type = "module",
                                            thresholds = seq(0, 0.02 , by = 0.001),
                                            export = TRUE,
                                            export_threshold = 0.068,
                                            prefix = "metaLup",
                                            file_format = ".tsv")
      
      # Downregulated Genes
      inmoduleDOWNmetaL <- metaL %>% filter(beta <= -1)
      genesDEGmetaLdown <- crossmap(inmoduleDOWNmetaL, id_col = "GeneSymbol", from = "symbol", to = "ensembl")
       # Mapped: 184 / 184  
      genesDEGmetaLdown <- intersect(rownames(genesDEGmetaLdown), rownames(TOM))
      length(genesDEGmetaLdown)
       # [1] 68
      cyt.metaLdown <- GenesetCytoscapeExport(geneset = genesDEGmetaLdown,
                                              TOM = TOM,
                                              blockwiseModules_network = bwnet,
                                              GeneAnnotation = GeneAnnot,
                                              geneset_type = "module",
                                              thresholds = seq(0, 0.02 , by = 0.001),
                                              export = TRUE,
                                              export_threshold = 0.010,
                                              prefix = "metaLdown",
                                              file_format = ".tsv")
      
    }
    
    # Exporting edge files for Non-Lesional vs Healthy DEGs
    metaNL <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\meta-analysis_NL_DEG.rds")
    {
      # Upregulated Genes
      inmoduleUPmetaNL <- metaNL %>% filter(beta >= 1)
      genesDEGmetaNLup <- crossmap(inmoduleUPmetaNL, id_col = "GeneSymbol", from = "entrez", to = "ensembl")
       # Mapped: 56 / 57
       # LOC105372047
      genesDEGmetaNLup <- intersect(rownames(genesDEGmetaNLup), rownames(TOM))
      length(genesDEGmetaNLup)
       # [1] 7
      cyt.metaNLup <- GenesetCytoscapeExport(geneset = genesDEGmetaNLup,
                                             TOM = TOM,
                                             blockwiseModules_network = bwnet,
                                             GeneAnnotation = GeneAnnot,
                                             geneset_type = "module",
                                             thresholds = seq(0, 0.1 , by = 0.001),
                                             export = TRUE,
                                             export_threshold = 0,
                                             prefix = "metaNLup",
                                             file_format = ".tsv")
      
      
      # Downregulated Genes
      inmoduleDOWNmetaNL <- metaNL %>% filter(beta <= -1)
      genesDEGmetaNLdown <- crossmap(inmoduleDOWNmetaNL, id_col = "GeneSymbol", from = "entrez", to = "ensembl")
       # Mapped: 12 / 13 
       # LOC112268102
      genesDEGmetaNLdown <- intersect(rownames(genesDEGmetaLdown), rownames(TOM))
      length(genesDEGmetaNLdown)
       # [1] 0
    }
      
    # Exporting edge files for Lesional Vs Non-Lesional
    metaLNL <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\meta-analysis_LNL_DEG_new.rds")
    {
      # Upregulated Genes
      inmoduleUPmetaLNL <- metaLNL %>% filter(beta >= 1)
      genesDEGmetaLNLup <- crossmap(inmoduleUPmetaLNL, id_col = "GeneSymbol", from = "symbol", to = "ensembl")
       # Mapped: 46 / 46 
      genesDEGmetaLNLup <- intersect(rownames(genesDEGmetaLNLup), rownames(TOM))
      length(genesDEGmetaLNLup)
       # [1] 29
      cyt.metaLNLup <- GenesetCytoscapeExport(geneset = genesDEGmetaLNLup,
                                             TOM = TOM,
                                             blockwiseModules_network = bwnet,
                                             GeneAnnotation = GeneAnnot,
                                             geneset_type = "module",
                                             thresholds = seq(0, 0.1 , by = 0.001),
                                             export = TRUE,
                                             export_threshold = 0.012,
                                             prefix = "metaLNLup",
                                             file_format = ".tsv")

        
      # Downregulated Genes
      inmoduleDOWNmetaLNL <- metaLNL %>% filter(beta <= -1)
      genesDEGmetaLNLdown <- crossmap(inmoduleDOWNmetaLNL, id_col = "GeneSymbol", from = "symbol", to = "ensembl")
       # Mapped: 413 / 465 
       # LOC100129034 PARD3-DT SLMO2-ATP5E LINC01698 LOC102723338 LOC102723803 LOC102724227 102724852 LOC105369579 
       # 105369903 LOC105370069 SLC35F4-AS1 LOC105370827 LOC105371988 LOC105372337 LOC105372934 LOC105373178 LOC105373563 
       # 105374013 LOC105374270 LOC105374294 LOC105374317 LOC105374491 LOC105375791 LOC105375925 LOC105375947 LOC105375951 
       # 105376074 LOC105376398 105376707 LOC105377847 LOC105377947 LOC105378685 LOC105378861 105378978 LOC107984005 107984204 
       # LOC107984452 107984774 LOC107985285 107986035 LOC107986058 107986170 LOC107986171 107986189 107986597 LOC107987112 
       # 112267876 112268014 112268284 H3P6 RPS14P3 
      genesDEGmetaLNLdown <- intersect(rownames(genesDEGmetaLNLdown), rownames(TOM))
      length(genesDEGmetaLNLdown)
       # [1] 193
      cyt.metaLNLdown <- GenesetCytoscapeExport(geneset = genesDEGmetaLNLdown,
                                              TOM = TOM,
                                              blockwiseModules_network = bwnet,
                                              GeneAnnotation = GeneAnnot,
                                              geneset_type = "module",
                                              thresholds = seq(0, 0.1 , by = 0.001),
                                              export = TRUE,
                                              export_threshold = 0.035,
                                              prefix = "metaLNLdown",
                                              file_format = ".tsv")
      
      
       }
    
  }



# 13. HUB Gene collection from Cytohubba
#--------------------------------------------------------------------------------
{
  # For module Lightcyan and midnightblue
  {
    # Locating module-specific CytoHubba hub-gene ranking files
    folder_path <- "D://DRYLAB//ATOPICDERMATITIS//Network Analysis//Cytoscape//Cytohubba Module Analysis"
    mod_folders <- list.dirs(folder_path, recursive = FALSE, full.names = TRUE)
    
    # Reataining only the Lightcyan and Midnightblue module directories
    mod_folders <- mod_folders[grepl(paste(c("Lightcyan", "Midnightblue"), collapse = "|"), basename(mod_folders))]
    hub_folders <- list.dirs(mod_folders, recursive = FALSE, full.names = TRUE)
    
    # Identifying the hub-gene ranking directories within each module
    hub_folders <- hub_folders[grepl("^Hub genes rank ?list$", basename(hub_folders), ignore.case = TRUE)]
    parent_modules <- basename(dirname(hub_folders))
    
    # Collecting CytoHubba hub-gene ranking lists
    module_hublist <- list()
    for (folder in hub_folders) {

      # Identify all CytoHubba ranking files
      rankfiles <- list.files(folder, pattern = "^Ranklist.*\\.csv$", full.names = TRUE)
      if (length(rankfiles) == 0) {
        cat("No ranklist found in:", folder, "\n")
        next
      }
      ranklist <- list()
      
      # Read individual CytoHubba ranking results
      for (file in rankfiles) {
        tryCatch({
          ranks <- read.csv(
            file,
            header = TRUE,
            fill = TRUE,
            check.names = FALSE,
            skip = 1
          )
          ranklist[[basename(file)]] <- ranks
        }, error = function(e) {
          cat("⚠️ Error reading:", file, "\n", e$message, "\n\n")
        })
      }
      
      # Store ranking lists according to their parent module
      module_name <- basename(dirname(folder))
      module_hublist[[module_name]] <- ranklist
    }
    
    
    # Identifying recurrent hub genes across CytoHubba ranking methods
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
    
    # Extracting common hub genes using a predefined recurrence threshold
    HUB_Genes <- Get_Common_HUB_Genes_flexible(module_hublist, cutoff = 8)
    
    # Converting module-specific hub-gene results to a unified table
    HUB_long <- map2_dfr( HUB_Genes,
                          names(HUB_Genes),
                          ~ mutate(.x, Module = .y)) %>%
                rename(HubGenes = Genes, Frequency = Freq)
    
    write.xlsx(HUB_long, "Common HUB Genes per module.xlsx", overwrite = TRUE)
    
  }
  
  # For meta-anlysis cytohubba genes
  {
    # Locating CytoHubba hub-gene ranking files
    folder_path_meta<- "D:\\DRYLAB\\ATOPICDERMATITIS\\meta-analysis\\Cytoscape meta"
    mod_folders_meta <- list.dirs(folder_path_meta, recursive = FALSE, full.names = TRUE)
    
    # Reataining only the contrast specific directories
    mod_folders_meta <- mod_folders_meta[grepl(paste(c("LNL", "L", "NL"), collapse = "|"), basename(mod_folders_meta))]
    hub_folders_meta <- list.dirs(mod_folders_meta, recursive = FALSE, full.names = TRUE)
    
    # Identifying the hub-gene ranking directories within each contrast specific directories 
    hub_folders_meta <- hub_folders_meta[ grepl("\\b(UP|DOWN)\\b", basename(hub_folders_meta), ignore.case = TRUE) ] 
    parent_modules_meta <- paste(basename(dirname(hub_folders_meta)), basename(hub_folders_meta))
    
    # Collecting CytoHubba hub-gene ranking lists
    module_hublist_meta <- list()
    for (folder in hub_folders_meta) {
      cat(folder, "\n")
      rankfiles_meta <- list.files(folder,
                                   pattern = ".*\\.csv$",
                                   full.names = TRUE)
      cat(rankfiles_meta)
      
      if (length(rankfiles_meta) == 0) {
        cat("No ranklist found in:", folder, "\n")
        next
      }
      
      ranklist <- list()
      
      for (file in rankfiles_meta) {
        
        tryCatch({
          
          ranks <- read.csv(
            file,
            header = TRUE,
            fill = TRUE,
            check.names = FALSE,
            skip = 1
          )
          
          ranklist[[basename(file)]] <- ranks
          
        }, error = function(e) {
          cat("⚠️ Error reading:", file, "\n", e$message, "\n\n")
        })
      }
      
      # ---- extract names ----
      cond <- basename(dirname(folder))   # L / LNL / NL
      cat(cond, "\n")
      direction   <- basename(folder)            # UP / DOWN
      cat(direction, "\n")
      
      # ---- nested assignment ----
      module_hublist_meta[[cond]][[direction]] <- ranklist
    }
    
    # Extracting common hub genes using a predefined recurrence threshold
    HUB_Genes_meta <- Get_Common_HUB_Genes_flexible(module_hublist_meta, cutoff = 8)
    
    # Converting module-specific hub-gene results to a unified table
    HUB_long_meta <- imap_dfr( HUB_Genes_meta,
                               ~ bind_rows(.x, .id = "Direction") %>%   # combine UP/DOWN
                                 mutate(Module = .y)) %>%
      rename(HubGenes = Genes, Frequency = Freq)
    
    write.xlsx(HUB_long_meta, "Common HUB Genes per Condition (meta).xlsx", overwrite = TRUE)
  }
}
  

 
# 12. HUB Gene identification
#--------------------------------------------------------------------------------
{
    # igraph methods : Degree + Betweenness + PowerCentrality + PageRank
    {
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
      # WGCNA modules
      modules_hub <- analyze_hub_modules(modules =  c("lightcyan", "midnightblue"),
                                                edge_dir = "D:\\DRYLAB\\ATOPICDERMATITIS\\Network Analysis",
                                                output_dir = "ModuleHubGenes",
                                                top_n = 10,
                                                methods = c("Degree", "Betweenness", "PowerCentrality", "PageRank"),
                                                hubcolor =c("#7AE2CF","#FF6969"), 
                                                nodecolor = "#F4CE14",
                                                plot_type = "ggraph",
                                                export_type = "png")
      
      # Meta-analysis
      igresmetaL <- analyze_hub_modules(modules =  c("metaLup", "metaLdown"),
                                                edge_dir = "D:\\DRYLAB\\ATOPICDERMATITIS\\Network Analysis",
                                                output_dir = "HubGeneResults(metaDEGsLnew)",
                                                top_n = 10,
                                                methods = c("Degree", "Betweenness", "PowerCentrality", "PageRank"),
                                                hubcolor =c("#7AE2CF","#FF6969"), 
                                                nodecolor = "#F4CE14",
                                                plot_type = "ggraph",
                                                export_type = "png")
      
      
      
      igresmetaNL <- analyze_hub_modules(modules =  c("metaNLup"),
                                             edge_dir = "D:\\DRYLAB\\ATOPICDERMATITIS\\Network Analysis",
                                             output_dir = "HubGeneResults(metaNLDEGs)",
                                             top_n = 10,
                                             methods = c("Degree", "Betweenness", "PowerCentrality", "PageRank"),
                                             hubcolor =c("#7AE2CF","#FF6969"), 
                                             nodecolor = "#F4CE14",
                                             plot_type = "ggraph",
                                             export_type = "png")
      
      
      
      igresmetaLNL <- analyze_hub_modules(modules =  c("metaLNLup", "metaLNLdown"),
                                              edge_dir = "D:\\DRYLAB\\ATOPICDERMATITIS\\Network Analysis",
                                              output_dir = "HubGeneResults(metaLNLDEGs)",
                                              top_n = 10,
                                              methods = c("Degree", "Betweenness", "PowerCentrality", "PageRank"),
                                              hubcolor =c("#7AE2CF","#FF6969"), 
                                              nodecolor = "#F4CE14",
                                              plot_type = "ggraph",
                                              export_type = "png")
      
      
         
        
      }
    
    
  
    # Intramodular connectivity (kWithin)
    {
      # Function for calculating (kWithin) for genesets or modules 
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
      # Hub gene from WGCNA modules
      modules_to_analyze <- c("lightcyan", "midnightblue")
      HUB_inConn <-  HUB_intramodularConnectivity( modules = modules_to_analyze,
                                                   module_membership_measure = module.membership.measure,
                                                   adjacency = adjacency,
                                                   wgcna_colors = bwnet$colors,
                                                   GeneSymbols = GeneAnnot )
      
      
      # Hub genes from DEG meta-analysis (meta-analysis) L
      metaDEGs_to_analyse_L <- list(genesDEGmetaLupnew, genesDEGmetaLdownnew)
      names(metaDEGs_to_analyse_L) <- c("genesDEGmetaLup", "genesDEGmetaLdown")
      HUB_inConnDEGsmetaL <-  HUB_intramodularConnectivity( modules = names(metaDEGs_to_analyse_L),
                                                                module_membership_measure = module.membership.measure,
                                                                adjacency = adjacency,
                                                                wgcna_colors = bwnet$colors,
                                                                GeneSymbols = GeneAnnot,
                                                                gene_sets = metaDEGs_to_analyse_L)
      
      
      # Hub genes from DEG meta-analysis (meta-analysis) NL
      metaDEGs_to_analyse_NL <-list(genesDEGmetaNLup, genesDEGmetaNLdown)
      names(metaDEGs_to_analyse_NL) <- c("genesDEGmetaNLup", "genesDEGmetaNLdown")
      HUB_inConnDEGsmetaNL <-  HUB_intramodularConnectivity( modules = names(metaDEGs_to_analyse_NL),
                                                             module_membership_measure = module.membership.measure,
                                                             adjacency = adjacency,
                                                             wgcna_colors = bwnet$colors,
                                                             GeneSymbols = GeneAnnot,
                                                             gene_sets = metaDEGs_to_analyse_NL)
      
      
      # Hub genes from DEG meta-analysis (meta-analysis) LNL
      metaDEGs_to_analyse_LNL <-list(genesDEGmetaLNLup, genesDEGmetaLNLdown)
      names(metaDEGs_to_analyse_LNL) <- c("genesDEGmetaLNLup", "genesDEGmetaLNLdown")
      HUB_inConnDEGsmetaLNL <-  HUB_intramodularConnectivity( modules = names(metaDEGs_to_analyse_LNL),
                                                              module_membership_measure = module.membership.measure,
                                                              adjacency = adjacency,
                                                              wgcna_colors = bwnet$colors,
                                                              GeneSymbols = GeneAnnot,
                                                              gene_sets = metaDEGs_to_analyse_LNL)
      
    } 
  }



