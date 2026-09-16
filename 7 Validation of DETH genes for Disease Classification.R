#===============================================================================
#         *** Validation of DETH genes for Disease Classification ***
#===============================================================================


 
#Packages
{
  library(dplyr)
  library(pROC)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(openxlsx) 
  library(viridis)
  library(uwot)
  library(ggpubr)
  library(ggnewscale)
  library(tidyverse)
  library(pheatmap)
  library(writexl)
  library(glmnet)
}



# Processing prior classification
#-------------------------------------------------------------------------------
{
  # Loading filtered, batch corrected integrated dataset and its metadata for validation
  FilteredDataset <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\FilteredDataset.rds")
  Metadata <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Validation\\Metadata.rds")
  
  # Mapping ensemble ids to gene symbols using predefined function crossmap
  FilteredDataset <- crossmap(FilteredDataset, from = "ensembl", to = "symbol")
  
  # Defining differentially expressed target hub (DETH) genes
  Tg <- c("PDZK1",	"INSIG1",	"FAR2",	"AADACL3",	"CTNNA1",	"PGK1",	"UBE2N",
          "PGAM1", "SERPINA12",	"CLDN1",	"FLG2",	"SCEL", "DLGAP5",	"NCAPG", 
          "CHI3L2",	"SPRR2E",	"LYPD2","OAS2",	"PARP9",	"OAS1",	"IFI44", 
          "MX1","CCL13",	"GALNT6",	"POLR3G",	"TTC39A",	"CCNB1",	"NELL2")
  
  # Extracting gene expression in lesional and healthy samples for DETH genes
  ExprTg <- as.matrix(FilteredDataset[intersect(rownames(FilteredDataset),Tg), 
                                      Metadata$Condition %in% c("Lesional", "Healthy")])
  Xtg <- t(ExprTg)
  
  # Define metadata for lesional/non-lesional versus healthy comparisons (Y)
  LesionalVsHealthy <- Metadata[ Metadata$Condition %in% c("Lesional", "Healthy"), ]
  Y <- factor(LesionalVsHealthy$Condition, levels = c("Healthy", "Lesional"))
  
  
  # Variables created after Elastic regressions and used for refined model
  genes.1se
  # lamda.1se
  # [1] "CHI3L2"    "TTC39A"    "OAS1"      "PGK1"      "NCAPG"     "POLR3G"   
  # [7] "CCNB1"     "SCEL"      "IFI44"     "GALNT6"    "FLG2"      "SERPINA12"
  # [13] "PDZK1"     "UBE2N"     "CCL13"     "AADACL3"   "LYPD2"     "SPRR2E" 
  Expr.1se <- as.matrix(FilteredDataset[genes.1se, Metadata$Condition %in% c("Lesional", "Healthy")])
  X.1se <- t(Expr.1se)
}



# 1. Heatmap 
#-------------------------------------------------------------------------------
{
  # Function to plot heat map for a given gene set from a dataset across 
  # skin conditions
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
  
  
  
  Target_hub_heat <- plotheatmap( FullExp = FilteredDataset,
                                  genes = Tg, 
                                  Metadata = Metadata, 
                                  clustering_method = "ward.D2",
                                  clustering_distance_rows = function(x) {
                                    as.dist(1 - cor(t(x), method = "pearson"))
                                  },
                                  clustering_distance_cols =function(x) {
                                    as.dist(1 - cor(t(x), method = "pearson"))
                                  },
                                  Plot_title = "Transcription Target Hub Gene Heatmap (Z-scored Expression)")
  
}



# 2. Logistic regression and ROC based evaluation
#-------------------------------------------------------------------------------
{
  # Univariate regression of DETH genes
  {
    # Estimating the association of each DETH genes with disease condition
    # using separate univariate logistic regression models
    uni_coeff_tg <- lapply(colnames(Xtg), function(g) {
      df <- data.frame( Y = Y,
                        X = Xtg[, g])
      fit <- glm(Y ~ X, data = df, family = binomial)
      coef(summary(fit))[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
    })
    names(uni_coeff_tg) <- colnames(Xtg)
    uni_coeff_tg <- do.call(rbind, uni_coeff_tg)
    uni_coeff_tg <- data.frame(Targets = colnames(Xtg), uni_coeff_tg)
    colnames(uni_coeff_tg) <- c( "Targets", "Estimate", "StdError", "Pvalue")
    uni_coeff_tg$CIlower <- uni_coeff_tg$Estimate - 1.96 * uni_coeff_tg$StdError
    uni_coeff_tg$CIupper <- uni_coeff_tg$Estimate + 1.96 * uni_coeff_tg$StdError
    uni_coeff_tg$OddsRatio <- exp(uni_coeff_tg$Estimate)
    write.xlsx(uni_coeff_tg, "Univariate Logistic Regression Coefficients of Targets new revised.xlsx", overwrite = TRUE)
    
    # Evaluating the discriminatory performance of individual genes using receiver 
    # operating characteristic (ROC) analysis and calculating the model-fit statistics
    uni_roc_list_TG <- lapply(colnames(Xtg), function(g) {
      
      # Predictor for this gene
      xg <- Xtg[, g]
      
      # Logistic regression for this gene
      fit <- glm(Y ~ xg, family = binomial)
      
      # McFadden R²
      LL_full <- as.numeric(logLik(fit))
      LL_null <- as.numeric(logLik(glm(Y ~ 1, family = binomial)))
      McFadden_R2 <- 1 - (LL_full / LL_null)
      
      CoxSnell_R2 <- 1 - exp((LL_null - LL_full) * 2 / n)
      Nagelkerke_R2 <- CoxSnell_R2 / (1 - exp(2 * LL_null / n))
      
      roc_obj <- roc(
        response = Y,
        predictor = xg,
        levels = c("Healthy","Lesional"))
      
      data.frame(Targets = g,
                 FPR = 1 - roc_obj$specificities,
                 TPR = roc_obj$sensitivities,
                 AUC = as.numeric(auc(roc_obj)),
                 R2 = McFadden_R2,
                 NR2 = Nagelkerke_R2,
                 Cxsnl = CoxSnell_R2)
    })
    uni_roc_list_TG <- bind_rows(uni_roc_list_TG)
    uni_auc_labels_TG <- uni_roc_list_TG %>%
      group_by( Targets ) %>% 
      summarize( AUC = unique(AUC),
                 R2  = unique(R2),
                 NR2 = unique(NR2),
                 CxsnlR = unique(Cxsnl))
    write.xlsx(uni_auc_labels_TG, "Univariate ROC-AUC Target new revised.xlsx", overwrite = TRUE)
    # Plotting the ROC curves for each DETH gene
    ggplot(uni_roc_list_TG, aes(x = FPR, y = TPR, color = Targets)) +
      geom_line(linewidth = 1.1, alpha = 0.85, color = "red") +
      geom_abline( linetype = "dashed", color = "black") +
      geom_text(  data = uni_auc_labels_TG,
                  aes( x = 0.3, y = 0.8,
                       label = paste0(  "AUC = ", round(AUC, 3), "\n",
                                        "R² = ", round(NR2, 3))),
                  color = "black", size = 6.3, inherit.aes = FALSE) +
      facet_wrap(~ Targets, scales = "free") +
      theme_bw(base_size = 14) +
      theme( panel.grid = element_blank(),
             strip.text = element_text(size = 18, face = "bold"),
             axis.text = element_text(color = "black", size = 16),
             axis.title = element_text(size = 24, face = "bold"),
             legend.position = "none",
             plot.title = element_text(size = 26, face = "bold", hjust = 0.5)) +
      labs( title = "Univariate ROC Curves for Targets",
            x = "False Positive Rate",
            y = "True Positive Rate")
    
  }
  
  # Multivariate regression of DETH genes 
  {
    # Evaluating the combined association of all DETH genes using a multivariate
    # logistic regression model
    LogiTg <- data.frame(Y, Xtg)
    multi_fit <- glm(Y ~ ., data = LogiTg, family = binomial)
    summary(multi_fit)
    multi_coeff <- as.data.frame(summary(multi_fit)$coefficients)
    multi_coeff$Targets <- rownames(multi_coeff)
    colnames(multi_coeff) <- c("Estimate", "StdError", "Z", "Pvalue", "Targets")
    multi_coeff <- multi_coeff[, c("Targets", "Estimate", "StdError", "Z", "Pvalue")]
    multi_coeff$CIlower <- multi_coeff$Estimate - 1.96 * multi_coeff$StdError
    multi_coeff$CIupper <- multi_coeff$Estimate + 1.96 * multi_coeff$StdError
    multi_coeff$OddsRatio <- exp(coef(multi_fit))
    write.xlsx( multi_coeff, "Multivariate Logistic Regression Coefficients of Targets new new.xlsx", overwrite = TRUE)
    multi_coeff <- multi_coeff[multi_coeff$Targets != "(Intercept)", ]
    
    
   # Point graph of Regression Coefficients with respective 95% confidence intervals
    ggplot(multi_coeff, aes(x = reorder(Targets, Estimate), y = Estimate)) +
      geom_point(aes(color = Targets), size = 3) +
      geom_errorbar(aes(ymin = CIlower, ymax = CIupper), width = 0.2, color = "steelblue") +
      scale_fill_manual(values = viridis(nrow(multi_coeff))) +
      geom_hline(yintercept = 0, linewidth = 0.8) +
      coord_flip() +
      theme_bw(base_size = 14) +
      labs(  title = "Logistic Regression Coefficients with 95% CI",
             y = "Effect Size (Log Odds)",
             x = "Gene") +
      theme( panel.grid = element_blank(),
             axis.text = element_text(color = "black"),
             axis.title = element_text(color = "black", face ="bold"),
             plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    
    
    
    # Evaluating the discriminatory performance of the DETH gene model
    # using receiver operating characteristic (ROC) analysis
    multi_pred <- predict(multi_fit, type = "response")
    multi_ROC <- roc(Y, multi_pred)
    auc(multi_ROC)                 # Area under the curve:0.9785
    # McFadden_R2
    LL_fullTG  <- as.numeric(logLik(multi_fit))
    LL_nullTG  <- as.numeric(logLik(glm(Y ~ 1, family = binomial)))
    multi_McFadden_R2_TG <- 1 - (LL_fullTG / LL_nullTG)
    multi_McFadden_R2_TG           # [1] 0.7376915
    # Nagelkerke R²
    n <- length(Y)
    Nagelkerke_R2_TG <- (1 - exp((LL_nullTG - LL_fullTG) * 2 / n)) / (1 - exp(LL_nullTG * 2 / n))
    Nagelkerke_R2_TG               # [1] 0.8531051
    multi_ROC_df <- data.frame( TPR = multi_ROC$sensitivities, FPR = 1 - multi_ROC$specificities)
    write_xlsx(multi_ROC_df, "multi-ROC Targets new revised.xlsx")
    
    
    # Plotting the ROC curves for the multivariate model with all 28 DETH genes
    ggplot(multi_ROC_df, aes(x = FPR, y = TPR)) +
      geom_line(color = "red", size = 1.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      annotate("text", x = 0.65, y = 0.2,
               label = paste0("AUC = ", round(auc(multi_ROC),4), "\n",
                              "McFadden R² = ", round(multi_McFadden_R2_TG,4), "\n",
                              "Nagelkerke R² = ", round(Nagelkerke_R2_TG, 4)),
               size = 5, hjust = 0, fontface = "bold") +
      theme_bw(base_size = 16) +
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs( title = "Multivariate Logistic Regression ROC",
            x = "False Positive Rate",
            y = "True Positive Rate")
    
  }
  
  # Multivariate Selected DETH genes at lambda.1se 
  {
    # Evaluating the combined association of 18 DETH genes using a multivariate
    # logistic regression model
    LogiTgselected <- data.frame(Y, X.1se)
    multi_fitL <- glm(Y ~ ., data = LogiTgselected, family = binomial)
    summary(multi_fitL)
    multi_coeffL <- as.data.frame(summary(multi_fitL)$coefficients)
    multi_coeffL$Targets <- rownames(multi_coeffL)
    colnames(multi_coeffL) <- c("Estimate", "StdError", "Z", "Pvalue", "Targets")
    multi_coeffL <- multi_coeffL[, c("Targets", "Estimate", "StdError", "Z", "Pvalue")]
    multi_coeffL$CIlower <- multi_coeffL$Estimate - 1.96 * multi_coeffL$StdError
    multi_coeffL$CIupper <- multi_coeffL$Estimate + 1.96 * multi_coeffL$StdError
    multi_coeffL$OddsRatio <- exp(coef(multi_fitL))
    write.xlsx( multi_coeffL, "Multivariate Logistic Regression Coefficients of Targets lambda.1se.xlsx", overwrite = TRUE)
    multi_coeffL <- multi_coeffL[multi_coeffL$Targets != "(Intercept)", ]
    
    # Point graph of Regression Coefficients with respective 95% confidence intervals
    ggplot(multi_coeffL, aes(x = reorder(Targets, Estimate), y = Estimate)) +
      geom_point(aes(color = Targets), size = 5) +
      geom_errorbar(aes(ymin = CIlower, ymax = CIupper), width = 0.2, color = "steelblue") +
      scale_fill_manual(values = viridis(nrow(multi_coeffL))) +
      geom_hline(yintercept = 0, linewidth = 0.8) +
      coord_flip() +
      theme_bw(base_size = 14) +
      labs(  title = bquote(bold("Lesional vs Healthy Logistic Regression Coefficients with 95% CI for Selected " ~lambda[1*SE] ~ " genes")),
             y = "Effect Size (Log Odds)",
             x = "Gene") +
      theme( panel.grid = element_blank(),
             axis.text = element_text(color = "black", size = 14),
             axis.title = element_text(color = "black", face ="bold", size = 18),
             plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
             legend.text = element_text(color = "black", size = 14, face = "bold"),
             legend.title = element_text(color = "black", size = 16, face = "bold") )
    
    
    # Evaluating the discriminatory performance of the DETH gene model
    # using receiver operating characteristic (ROC) analysis
    multi_predL <- predict(multi_fitL, type = "response")
    multi_ROCL <- roc(Y, multi_predL)
    auc(multi_ROCL)                             # Area under the curve: 0.9764
    # McFadden_R2
    LL_fullL  <- as.numeric(logLik(multi_fitL))
    LL_nullL <- as.numeric(logLik(glm(Y ~ 1, family = binomial)))
    multi_McFadden_R2_L <- 1 - (LL_fullL / LL_nullL)
    multi_McFadden_R2_L                         # [1] 0.7218952
    # Nagelkerke R²
    n <- length(Y)
    Nagelkerke_R2_L <- (1 - exp((LL_nullL - LL_fullL) * 2 / n)) / (1 - exp(LL_nullL * 2 / n))
    Nagelkerke_R2_L                             # [1] 0.8424514
    multi_ROCL_df<- data.frame( TPR = multi_ROCL$sensitivities, 
                                FPR = 1 - multi_ROCL$specificities)
    write_xlsx(multi_ROCL_df, "multi-ROC_Lambda.min Targets.1se.xlsx")
    
    
    # Plotting the ROC curves for the multivariate model with 18 DETH genes
    ggplot(multi_ROCL_df, aes(x = FPR, y = TPR)) +
      geom_line(color = "red", size = 1.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      annotate("text", x = 0.65, y = 0.2,
               label = paste0("AUC = ", round(auc(multi_ROCL),4), "\n",
                              "McFadden R² = ", round(multi_McFadden_R2_L,4), "\n",
                              "Nagelkerke R² = ", round(Nagelkerke_R2_L, 4)),
               size = 5, hjust = 0, fontface = "bold") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 16, face ="bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
      labs( title = bquote(bold("Multivariate Logistic Regression ROC for Selected"~lambda[1*SE] ~ " genes")),
            x = "False Positive Rate",
            y = "True Positive Rate")
    
  }

}
 


# 3. Elastic Net regression to obtain a refined model
#-------------------------------------------------------------------------------
{
  # Fitting a cross-validated Elastic Net logistic regression model to identify
  # genes associated with disease status using binomial family of regerssion for
  # binary disease classification
  Elastic_tg <- cv.glmnet( x = Xtg,
                           y = Y,
                           family = "binomial",   
                           alpha = 0.5  )            
  # Visualising the cross-validation error across the regularization path
  plot(Elastic_tg)
  # Number of non-zero coefficients at the two selected regularization parameters
  Elastic_tg$nzero[Elastic_tg$index[1]]  # lambda.min s48 17
  Elastic_tg$nzero[Elastic_tg$index[2]]  # lambda.1se s37 18
  
  
  # Extracting model coefficients at lambda.1se
  Elastic_tg_coef_list1.se <- coef(Elastic_tg, s = "lambda.1se")
  Elastic_tg_coef_list1.se <- as.matrix(Elastic_tg_coef_list1.se)
  Elastic_tg_coef_list1.se <- data.frame( Gene = rownames(as.matrix(Elastic_tg_coef_list1.se)),
                                          Coefficient = as.numeric(Elastic_tg_coef_list1.se))
  rownames(Elastic_tg_coef_list1.se) <- Elastic_tg_coef_list1.se$Gene
  write_xlsx(Elastic_tg_coef_list1.se, "Elastic_Regression_Coefficients.1se.xlsx")
  genes.1se <- rownames(Elastic_tg_coef_list1.se)[Elastic_tg_coef_list1.se$Coefficient != 0]
  genes.1se <- setdiff(genes.1se, "(Intercept)")
  genes.1se
  # [1] "CHI3L2"    "TTC39A"    "OAS1"      "PGK1"      "NCAPG"     "POLR3G"   
  # [7] "CCNB1"     "SCEL"      "IFI44"     "GALNT6"    "FLG2"      "SERPINA12"
  # [13] "PDZK1"     "UBE2N"     "CCL13"     "AADACL3"   "LYPD2"     "SPRR2E" 
  
  
  # Extracting model coefficients at lambda.min
  Elastic_tg_coef_list.min <- coef(Elastic_tg, s = "lambda.min")
  Elastic_tg_coef_list.min <- as.matrix(Elastic_tg_coef_list.min)
  Elastic_tg_coef_list.min <- data.frame( Gene = rownames(as.matrix(Elastic_tg_coef_list.min)),
                                          Coefficient = as.numeric(Elastic_tg_coef_list.min))
  rownames(Elastic_tg_coef_list.min) <- Elastic_tg_coef_list.min$Gene
  write_xlsx(Elastic_tg_coef_list.min, "Elastic_Regression_Coefficients.min_new revised.xlsx")
  genes.min <- rownames(Elastic_tg_coef_list.min)[Elastic_tg_coef_list.min$Coefficient != 0]
  genes.min <- setdiff(genes.min, "(Intercept)")
  genes.min
  # [1] "CTNNA1"    "CHI3L2"    "TTC39A"    "PGK1"      "NCAPG"     "POLR3G"   
  # [7] "CCNB1"     "SCEL"      "IFI44"     "FLG2"      "SERPINA12" "PDZK1"    
  # [13] "UBE2N"     "CCL13"     "AADACL3"   "LYPD2"     "SPRR2E"   
  
}
  


# 4. Principel Component Analysis
#-------------------------------------------------------------------------------
{
  # Function to perform principle component analysis and visualising the results
  # for a given gene set from a dataset across skin conditions
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
  
  
  # PCA of 28 DETH genes
  pca_tg_rev <- pca_val(Expression = ExprTg,
                                Metadata = Metadata, 
                                colour_guide = NULL,
                                PCs = c("PC1", "PC2"),
                                title = bquote(bold("PCA of 28 Differentially expressed Target HUB genes")))
  
  # PCA of 18 DETH genes selected by elastic net lamda 1.se
  pca_tg_1se <- pca_val(Expression = Expr.1se,
                                Metadata = Metadata, 
                                colour_guide = NULL,
                                PCs = c("PC1", "PC2"),
                                title = bquote(bold("PCA plot of Selected " ~lambda[1*SE] ~ " Differentially expressed Target HUB genes")))
  
    }



# 5. Uniform Manifold Projection Approximation
#-------------------------------------------------------------------------------
{
  
  # Function to perform uniform manifold approximation projection and 
  # visualising the resultsfor a given gene set from a dataset across skin 
  # conditions
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
  
  # UMAP of 28 DETH genes
  UMAPtg <- plot_umap(expression_data = ExprTg ,
                                  colour = c( "Healthy" = "#4f6367", "Lesional" = "#0077c8"),
                                  Metadata = Metadata,
                                  title = bquote(bold("UMAP of 28 Differentially expressed Target HUB genes")))
  
  # UMAP of 18 DETH genes selected by elastic net lamda 1.se
  UMAPtg1.se <- plot_umap(expression_data = Expr.1se ,
                                  colour = c( "Healthy" = "#4f6367", "Lesional" = "#0077c8"),
                                  Metadata = Metadata,
                                  title = bquote(bold("UMAP of Selected " ~lambda[1*SE] ~ " Differentially expressed Target HUB genes")))
    
  }
  

