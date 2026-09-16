#===============================================================================
#                *** Differentially Regulated Link Analysis ***            
#-------------------------------------------------------------------------------
#                                *** GENIE3 ***                     
#===============================================================================


#Packages
{
  library(GENIE3)
  library(dorothea)
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(doRNG)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(writexl)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
}


# 1. Data Preprocessing
#-------------------------------------------------------------------------------
{
  # Loading the batch corrected master expression dataframe and transpose it to 
  # obtain samples as rows and genes as columnsa as the input orientation required for WGCNA.
  MasterAD <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Network Analysis\\Master_Dataset.rds")
  tMasterAD <- t(MasterAD)
  dim(tMasterAD)
  # [1]   826 13859
  
  # Normalising the gene expression to scale
  TsclExpr <- scale(tMasterAD)
  sclExpr <- t(TsclExpr)
  saveRDS(sclExpr, "Scaled Master_Dataset.rds")
  
  # Loading the metadata
  metadata <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\Network Analysis\\masterMetadata.rds")
  
  # Gene Annotation
  GeneSymbols <- readRDS("D:\\DRYLAB\\ATOPICDERMATITIS\\GeneSymbols.rds")
  
}



# 2: Filtering genes with low expression variance (bottom 25%)
#-------------------------------------------------------------------------------
{
   # Calculating variance 
  expr_var <- apply(sclExpr, 1, var)
  
  # 25th quantile
  variance_threshold <- quantile(expr_var, 0.25)
  
  # Retaining genes with variance greater than 25th quantile 
  sclExpr_fil <- sclExpr[expr_var > variance_threshold, ]
  dim(sclExpr_fil)
  # [1] 9322  826
}



# 3: Defining transcription factor list from DoRothEA Database 
#-------------------------------------------------------------------------------
{
  # Retrieving the human transcription factor–target regulatory information from DoRothEA.
  dorothea_data <- dorothea::dorothea_hs
  
  # Unique transcription factor gene symbols
  tf_symbols <- unique(dorothea_data$tf)
  length(tf_symbols)
  # [1] 1333
  
  # TFs that are present in gene annotation data
  tf_info <- GeneSymbols[GeneSymbols$hgnc_symbol %in% tf_symbols, c("ensembl_gene_id", "hgnc_symbol")]
  nrow(tf_info)
  # [1] 1330
  
  # TFs that are present in gene expression data
  sum(tf_info$ensembl_gene_id %in% rownames(sclExpr_fil))
  # [1] 772
  
  # Retaining only TFs found in expression data
  tf_table <- tf_info[tf_info$ensembl_gene_id %in% rownames(sclExpr_fil), ]
  
}



# 4: GENIE3 network in skin conditions
#-------------------------------------------------------------------------------
{
  # Defining samples and expression for each skin condition
  Healthy <- rownames(metadata)[metadata$Condition == "Healthy"]
  exprH <- sclExpr_fil[, Healthy]
  
  Lesional <- rownames(metadata)[metadata$Condition == "Lesional"]
  exprL <- sclExpr_fil[, Lesional]
  
  Nonlesional <- rownames(metadata)[metadata$Condition == "Non-lesional"]
  exprNL <- sclExpr_fil[, Nonlesional]
  
  
  # Performing GENIE3 for each skin condition with predefined TF table and fixed
  # seed to ensure reproducibility
  set.seed(123)
  netH <- GENIE3(exprH, regulators = tf_table$ensembl_gene_id, nCores = 8, verbose = TRUE)
  
  set.seed(123)
  netL <- GENIE3(exprL, regulators = tf_table$ensembl_gene_id, nCores = 8, verbose = TRUE)
  
  set.seed(123)
  netNL <- GENIE3(exprNL, regulators = tf_table$ensembl_gene_id, nCores = 8, verbose = TRUE)
  
  
  # Collecting the results of GENIE3 in a list
  GRNs <- list(Healthy = netH,
               Lesional = netL,
               Nonlesional = netNL)
  
  
  # Processing the GENIE3 results to retain only high confidence links in each skin condition
  Cond.Links <- lapply(GRNs, function(x){
    df <- x 
    links <- getLinkList(df)
    
    # Mapping the ensemble identifiers to gene symbols
    links$regulatorSymbol <- GeneSymbols$hgnc_symbol[match(links$regulatoryGene, GeneSymbols$ensembl_gene_id)]
    links$targetSymbol <- GeneSymbols$hgnc_symbol[match(links$targetGene, GeneSymbols$ensembl_gene_id)]
    
    # Determining threshold for high confidence i.e., top 1% of interaction weights
    weight_threshold <- quantile(links$weight, 0.99)
    weight_threshold
    high_conf <- links %>% filter(weight >= weight_threshold)
    
    # Standardised representation of each link with 6 only attributes
    high_conf$RegulatoryPair <- paste(high_conf$regulatorSymbol, high_conf$targetSymbol, sep= " ==> " )
    high_conf <- high_conf[, c("regulatoryGene", "regulatorSymbol", 
                               "targetGene", "targetSymbol", "weight", "RegulatoryPair" )]
    
    # Calculating the total regulatory influence of a TF as summary 
    tf_sum <- high_conf %>%
      group_by(regulatorSymbol) %>%
      summarise(sum_weight = sum(weight),
                n_targets = n(),
                .groups = 'drop')
    
    
    return(list(tf_summary = tf_sum,
                high_conf_links = high_conf))
  })
  
  
  ## Performing GENIE3 overall conditions with predefined TF table and fixed
  GENnet <- GENIE3(sclExpr_fil, regulators = tf_table$ensembl_gene_id, nTrees = 1000, K = "sqrt", nCores = 8, verbose = TRUE) 
  Links <- getLinkList(GENnet)
  Links$regulatorSymbol <- GeneSymbols$hgnc_symbol[  match(Links$regulatoryGene, GeneSymbols$ensembl_gene_id)]
  Links$targetSymbol <- GeneSymbols$hgnc_symbol[match(Links$targetGene, GeneSymbols$ensembl_gene_id)]
  Links <- Links[, c("regulatoryGene", "regulatorSymbol", "targetGene", "targetSymbol", "weight")]
  Links$RegulatoryPair <- paste(Links$regulatorSymbol," ==> ", Links$targetSymbol)
  weight_threshold <- quantile(Links$weight, 0.99)
  Highlinks <- Links %>% filter(weight >= weight_threshold)
  
}



# 5: Transcription factor summary
#-------------------------------------------------------------------------------
{
  # Appending condition labels to the TF summary results to facilitate integration and comparison across the three skin conditions
  colnames(Cond.Links$Healthy$tf_summary)[-1]  <- paste0(colnames(Cond.Links$Healthy$tf_summary)[-1], "_Healthy")
  colnames(Cond.Links$Lesional$tf_summary)[-1] <- paste0(colnames(Cond.Links$Lesional$tf_summary)[-1], "_Lesional")
  colnames(Cond.Links$Nonlesional$tf_summary)[-1] <- paste0(colnames(Cond.Links$Nonlesional$tf_summary)[-1], "_Nonlesional")
  
  # Integrating  transcription factor regulatory influence across skin conditions
  Cond.comparison <- merge( Cond.Links$Healthy$tf_summary, Cond.Links$Lesional$tf_summary, by = "regulatorSymbol", all = TRUE)
  Cond.comparison <- merge( Cond.comparison, Cond.Links$Nonlesional$tf_summary, by = "regulatorSymbol", all = TRUE)
  Cond.comparison[is.na(Cond.comparison)] <- 0
  Cond.long <- Cond.comparison %>%
               pivot_longer(cols = -regulatorSymbol,
                            names_to = c(".value", "Condition"),
                            names_pattern = "(.+)_(Healthy|Lesional|Nonlesional)")
  
  # Quantify differential transcription factor regulatory influence
  # magnitude and direction of difference
  Cond.comparison$wdiff_LvH <- Cond.comparison$sum_weight_Lesional - Cond.comparison$sum_weight_Healthy
  Cond.comparison$wdiff_NLvH <- Cond.comparison$sum_weight_Nonlesional - Cond.comparison$sum_weight_Healthy
  Cond.comparison$wdiff_LvNL <- Cond.comparison$sum_weight_Lesional - Cond.comparison$sum_weight_Nonlesional
  # log fold change of transcription factor weight
  Cond.comparison$wLFC_LvH <- log2((Cond.comparison$sum_weight_Lesional + 1) / (Cond.comparison$sum_weight_Healthy + 1))
  Cond.comparison$wLFC_NLvH <- log2((Cond.comparison$sum_weight_Nonlesional + 1) / (Cond.comparison$sum_weight_Healthy + 1))
  Cond.comparison$wLFC_LvNL <- log2((Cond.comparison$sum_weight_Lesional + 1) / (Cond.comparison$sum_weight_Nonlesional + 1))
  
  write_xlsx(Cond.comparison, "Transcription Factor Condition-Comparison.xlsx")
 
  # Pairwise Scatter plot of regulatory influence of TFs
  # Lesinoal Vs Healthy
  {
    ggplot(Cond.comparison,
           aes(x = sum_weight_Healthy,
               y = sum_weight_Lesional,
               color = wLFC_LvH)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = regulatorSymbol),
                      size = 4,
                      color = "black") +
      scale_color_gradient2(
        low = "#0b5351",
        mid = "#fed766",
        high = "#d62828",
        midpoint = 0,
        name = "Log2 FoldChange") +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray") +
      theme_classic() +
      labs(
        title = "Transcription Factor Regulatory Influence",
        subtitle = "Comparison between Lesional and Healthy",
        x = "Total Regulatory Weight (Healthy)",
        y = "Total Regulatory Weight (Lesional)"
      ) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size =16),
            axis.title = element_text (face = "bold", size = 12),
            plot.subtitle = element_text(face = "italic", hjust = 0.5, size =14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))
    
  }
  
  # Non-Lesional Vs Healthy
  {
    ggplot(Cond.comparison,
           aes(x = sum_weight_Healthy,
               y = sum_weight_Nonlesional,
               color = wLFC_NLvH)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = regulatorSymbol),
                      size = 4,
                      color = "black") +
      scale_color_gradient2(
        low = "#0b5351",
        mid = "#fed766",
        high = "#d62828",
        midpoint = 0,
        name = "Log2 FoldChange") +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray") +
      theme_classic() +
      labs(
        title = "Transcription Factor Regulatory Influence",
        subtitle = "Comparison between Non-Lesional and Healthy",
        x = "Total Regulatory Weight (Healthy)",
        y = "Total Regulatory Weight (Non-Lesional)"
      ) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size =16),
            axis.title = element_text (face = "bold", size = 12),
            plot.subtitle = element_text(face = "italic", hjust = 0.5, size =14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))
    
  }
  
  #Lesional Vs Non-Lesional
  {
    ggplot(Cond.comparison,
           aes(x = sum_weight_Nonlesional,
               y = sum_weight_Lesional,
               color = wLFC_LvNL)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = regulatorSymbol),
                      size = 4,
                      color = "black") +
      scale_color_gradient2( low = "#0b5351",
                             mid = "#fed766",
                             high = "#d62828",
                             midpoint = 0,
                             name = "Log2 FoldChange") +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray") +
      theme_classic() +
      labs( title = "Transcription Factor Regulatory Influence",
            subtitle = "Comparison between Lesional and Non-Lesional",
            x = "Total Regulatory Weight (Non-Lesional)",
            y = "Total Regulatory Weight (Lesional)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size =16),
            axis.title = element_text (face = "bold", size = 12),
            plot.subtitle = element_text(face = "italic", hjust = 0.5, size =14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))
    
  }
  
  # Heat map of LFC_sum_weightage of TFs
  {
    # Heatmap of top 10 genes from all contrasts
    
    top10TFs <- c(  "FOXM1",	"ARID5A",	"BBX",	"BHLHE40",	"CDC5L",	"CEBPZ",	"CREBL2",	"CREM",	"DMBX1",	"E2F7",
                    "EGR2",	"ETV1",	"GATA2",	"GTF2B",	"HHEX",	"HINFP",	"IRF4",	"IRF6",	"IRX1",	"IRX6", "MYBL2",
                    "NR1H3",	"PPARG",	"RARG",	"RBPJ",	"SNAI2",	"SNAI3",	"SP3",	"TCF12",	"THAP11", "ZHX1",
                    "ZNF117",	"ZNF254",	"ZNF358",	"ZNF511",	"ZNF561",	"ZNF574",	"ZNF580",	"ZNF584",	"ZNF672", "KLF10",
                    "MAX", "MSC",	"ZBTB7C")
    top10LFC <- Cond.comparison[Cond.comparison$regulatorSymbol %in% top10TFs, c(1,2,4,6)]
    rownames(top10LFC) <- top10LFC$regulatorSymbol
    top10LFC <- top10LFC %>% rename( Healthy = sum_weight_Healthy,
                                     `Non-Lesional` = sum_weight_Nonlesional,
                                     Lesional = sum_weight_Lesional)
    top10LFC <- top10LFC[ , c(2,4,3)]
    top10LFC <- as.matrix(top10LFC)
    
    # Removing rows with no value at all
    top10LFC <- top10LFC[rowSums(abs(top10LFC)) > 0, ]
    
    # Scaling
    top10LFC <- t(scale(t(top10LFC)))
    
    # Colors
    range(top10LFC)
    # [1] -1.151451  1.154680
    col_fun <- colorRamp2(
      c(-0.2060594 , 20.8315811),
      c( "white", "#B2182B" ))
    
    # Export PNG
    tiff( "sumWeightage_Heatmap no clustering.tiff", width = 5.1, height = 8.50, units = "in", res = 1200)
    
    
    ht <- Heatmap(
      top10LFC,
      col = col_fun,
      
      # clustering
      cluster_rows =TRUE,
      cluster_columns = FALSE,
      
      # row labels
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 10, fontface = "bold"),
      
      # ---- COLUMN LABELS CENTERED ----
      column_names_side = "bottom",
      column_names_centered = TRUE,
      column_names_rot = 0,
      column_names_gp = gpar(fontsize = 10, fontface = "bold"),
      
      
      rect_gp = gpar(col = NA),
      
      name = "Sum Weightage",
      
      # legend styling
      heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"),
                                  labels_gp = gpar(fontsize = 9)),
      
      
      row_title = "Transcription Factors",
      row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    )
    
    draw(
      ht,
      column_title = "Sum Weightage of Transcription Factors",
      column_title_gp = gpar(fontsize = 12, fontface = "bold")
    )
    #ht
    dev.off()
  }
}



# 6: Skin condition-specific links
#-------------------------------------------------------------------------------
{
    # Regulatory pairs exclusively found only in specific skin conditions
    H_pair <- setdiff(Cond.Links$Healthy$high_conf_links$RegulatoryPair, union(Cond.Links$Lesional$high_conf_links$RegulatoryPair, Cond.Links$Nonlesional$high_conf_links$RegulatoryPair))
    H_pair <- Cond.Links$Healthy$high_conf_links[Cond.Links$Healthy$high_conf_links$RegulatoryPair %in% H_pair, ]
    
    L_pair <- setdiff(Cond.Links$Lesional$high_conf_links$RegulatoryPair, union(Cond.Links$Nonlesional$high_conf_links$RegulatoryPair, Cond.Links$Healthy$high_conf_links$RegulatoryPair))
    L_pair <- Cond.Links$Lesional$high_conf_links[Cond.Links$Lesional$high_conf_links$RegulatoryPair %in% L_pair, ]
    
    NL_pair <- Cond.Links$Nonlesional$high_conf_links[Cond.Links$Nonlesional$high_conf_links$RegulatoryPair %in% NL_pair, ]
    NL_pair <- setdiff(Cond.Links$Nonlesional$high_conf_links$RegulatoryPair, union(Cond.Links$Lesional$high_conf_links$RegulatoryPair, Cond.Links$Healthy$high_conf_links$RegulatoryPair))

    
    # Retrieving symbol annotations for those missing
    syngo <- read_excel("C:\\Users\\HP\\Downloads\\syngo_id_convert_2026-02-24_12-18.xlsx")
    H_pair <- H_pair %>% left_join( syngo %>% dplyr::select(query, symbol), by = c("targetGene" = "query")) %>%
                         mutate(targetSymbol = coalesce(symbol, targetSymbol)) %>%
                         dplyr::select(-symbol) %>%
                         mutate(RegulatoryPair = paste(regulatorSymbol, targetSymbol, sep= " ==> " ))
    L_pair <- L_pair %>% left_join( syngo %>% dplyr::select(query, symbol), by = c("targetGene" = "query")) %>%
                        mutate(targetSymbol = coalesce(symbol, targetSymbol)) %>%
                        dplyr::select(-symbol) %>%
                        mutate(RegulatoryPair = paste(regulatorSymbol, targetSymbol, sep= " ==> " ))
    NL_pair <- NL_pair %>% left_join( syngo %>% dplyr::select(query, symbol), by = c("targetGene" = "query")) %>%
                           mutate(targetSymbol = coalesce(symbol, targetSymbol)) %>%
                           dplyr::select(-symbol) %>%
                          mutate(RegulatoryPair = paste(regulatorSymbol, targetSymbol, sep= " ==> " ))
    
    
    # Exporting the condition-specific targets, TFs and links
    UniTarget <- list(  HealthyTg    = data.frame(targetSymbol = unique(H_pair$targetSymbol)),
                        LesionalTg   = data.frame(targetSymbol = unique(L_pair$targetSymbol)),
                        NonLesionalTg = data.frame(targetSymbol = unique(NL_pair$targetSymbol)))
    write_xlsx(UniTarget, "Condition-Wise-Targets.xlsx")
    
    
    UniTfs <- list( HealthyTF    = data.frame(regulatorySymbol = unique(H_pair$regulatorSymbol)),
                    LesionalTF  = data.frame(regulatorySymbol = unique(L_pair$regulatorSymbol)),
                    NonLesionalTF = data.frame(regulatorySymbol = unique(NL_pair$regulatorSymbol)))
    write_xlsx(UniTfs, "Condition-Wise-TFs.xlsx")
    
    
    UniTftgs <- list( HealthyTFtg    = data.frame(RegulatoryPair = unique(H_pair$RegulatoryPair)),
                      LesionalTFtg  = data.frame(RegulatoryPair = unique(L_pair$RegulatoryPair)),
                      NonLesionalTtgF = data.frame(RegulatoryPair = unique(NL_pair$RegulatoryPair)))
    write_xlsx(UniTftgs, "Condition-Wise-TFtgs.xlsx")
    
  }



# 7: Network visualization
#-------------------------------------------------------------------------------
{
  # Creating a network graph from top links
  # Use top 500 links for visualization
  top_H_links <- Cond.Links$Healthy$high_conf_links[1:500, ]
  top_L_links <- Cond.Links$Lesional$high_conf_links[1:500, ]
  top_NL_links <- Cond.Links$Nonlesional$high_conf_links[1:500, ]
  
  # Compiling the high confidence links into list
  linklist <- list(H = top_H_links,
                   NL = top_NL_links,
                   L = top_L_links)
  
  # igraph object
  H_g <- graph_from_data_frame( top_H_links[, c("regulatorSymbol", "targetSymbol", "weight")], directed = TRUE)
  NL_g <- graph_from_data_frame( top_NL_links[, c("regulatorSymbol", "targetSymbol", "weight")], directed = TRUE)
  L_g <- graph_from_data_frame( top_L_links[, c("regulatorSymbol", "targetSymbol", "weight")], directed = TRUE)
  
  # Compiling the igraph  object into list
  Iobj <- list(H = H_g,
               NL= NL_g,
               L = L_g)
  
  # Network attributes and plotting
  Networks <- mapply(function(g, links, name){
    
    # ---- Node properties ----
    V(g)$degree     <- degree(g, mode = "all")
    V(g)$in_degree  <- degree(g, mode = "in")
    V(g)$out_degree <- degree(g, mode = "out")
    
    # ---- TF properties ----
    is_tf <- V(g)$name %in% unique(links$regulatorSymbol)
    V(g)$type <- ifelse(is_tf, "TF", "Target")
    
    # ---- Visual attributes ----
    V(g)$color <- ifelse(V(g)$type == "TF", "#07553B", "#CED46A")
    V(g)$size  <- sqrt(V(g)$degree) * 3 + 3
    
    # ---- Edge attributes ----
    E(g)$width <- (E(g)$weight / max(E(g)$weight)) * 3
    
    # ---- Saving plot ----
    png(paste0(name, ".png"), width = 1200, height = 1200, res = 100)
    par(mar = c(1, 1, 3, 1))
    set.seed(123)
    
    plot <- plot(g,
         vertex.label = ifelse( V(g)$degree > quantile(V(g)$degree, 0.9), V(g)$name, NA),
         vertex.label.cex = 0.7,
         vertex.label.color = ifelse( V(g)$type == "TF", "white", "#07553B"),,
         edge.arrow.size = 0.3,
         edge.color = "black",
         layout = layout_with_fr(g),
         main = paste("Gene Regulatory Network -", name))
    legend("topright",
           legend = c("Transcription Factor", "Target Gene"),
           col = c("#07553B", "#CED46A"),
           pch = 19,
           cex = 0.8,
           bty = "n")
    dev.off()
 
    # Save network statistics
    network_stats <- data.frame(Metric = c( "Nodes", "Edges", "Density", "TFs", "Targets" ),
                                Value = c( vcount(g),
                                           ecount(g),
                                           edge_density(g),
                                           sum(V(g)$type == "TF"),
                                           sum(V(g)$type == "Target")))
    
    return(list(graph = g,
                stats = network_stats))
     
    }, Iobj, linklist, names(Iobj), SIMPLIFY = FALSE)
  names(Networks) <- names(Iobj)

  # Exporting network stats in excel
  lapply(names(Networks), function(name) {
    network_stats <- Networks[[name]]$stats
    writexl::write_xlsx(network_stats,paste0("network_stats_", name, ".xlsx")) })
  
 
}



# 8: Sub Network analysis
#-------------------------------------------------------------------------------
{
  # Feed-forward loops: TF1 → TF2 → Target, TF1 → Target
  # These represent coordinated regulation
  # Subnetwork Attributes and Plotting
  Subnetworks <- lapply(names(Networks), function(name){
    
    g <- Networks[[name]]$graph
    
    # ---- TF-only subgraph ----
    tf_nodes <- V(g)[type == "TF"]$name
    tf_subgraph <- igraph::induced_subgraph(g, tf_nodes)
    
    # Removing isolated TFs
    tf_subgraph <- delete_vertices(tf_subgraph, V(tf_subgraph)[igraph::degree(tf_subgraph) == 0])
    
    # ---- Subnetwork ----
    if (ecount(tf_subgraph) > 0) {
      
      # out-degree
      V(tf_subgraph)$out_degree <- igraph::degree(tf_subgraph, mode = "out")
      
      # Keeping top 50% TFs 
      threshold <- quantile(V(tf_subgraph)$out_degree, 0.0)
      keep_nodes <- V(tf_subgraph)[ out_degree >= threshold]$name
      tf_subgraph <- igraph::induced_subgraph(tf_subgraph, keep_nodes)
      
      # Collapsing strongly connected components
      comp <- igraph::components(tf_subgraph, mode = "strong")
      
      # Contracting graph
      condensed <- igraph::contract(tf_subgraph, comp$membership)
      
      # Removing multi-edges
      condensed <- igraph::simplify(condensed)
      
      # hierarchy depth on condensed DAG
      roots <- V(condensed)[igraph::degree(condensed, mode = "in") == 0]
      dist_matrix <- distances(condensed, v = roots, mode = "out")
      depth_meta <- apply(dist_matrix, 2, function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0) return(0)
        min(x)
      })
      
      # Mapping meta-node depth back to original nodes
      node_depth <- depth_meta[comp$membership]
      V(tf_subgraph)$depth <- node_depth
      
      # Sugiyama layout (hierarchical)
      sug <- layout_with_sugiyama(tf_subgraph, hgap = 20, vgap = 20, maxiter = 500,)
      lay <- sug$layout
      layers <- split(seq_len(nrow(lay)),round(lay[,2], 3))
      for(idx in layers){
        if(length(idx) > 1){
          
          ord <- order(lay[idx,1])
          spacing <- 2000
          lay[idx[ord],1] <- seq( -(length(idx)-1)/10*spacing, 
                                  (length(idx)-1)/10*spacing,
                                  length.out = length(idx))
          centre <- mean(lay[idx,1])
          lay[idx,1] <- centre + (lay[idx,1] - centre) * 1.5
        }
      }
      lay[,1] <- lay[,1] * 2.5
      
      # Layers based on depth ranking
      ranked_depth <- rank(V(tf_subgraph)$depth, ties.method = "first")
      V(tf_subgraph)$layer <- cut( ranked_depth, breaks = 3, labels = c("Layer 1", "Layer 2", "Layer 3"))
      
      # Dynamic layer coloring
      layer_levels <- levels(factor(V(tf_subgraph)$layer))
      n_layers <- length(layer_levels)
      
      # Generate colors dynamically
      base_colors <- c( "#EDCD44","#B6E696", "#CADCFC")
      layer_colors <- setNames(base_colors[seq_len(min(n_layers, length(base_colors)))],
                               layer_levels)
      V(tf_subgraph)$color <- layer_colors[as.character(V(tf_subgraph)$layer)]
      
      # Vertex size
      lab_len <- nchar(V(tf_subgraph)$name)
      V(tf_subgraph)$size <-   pmax(16, lab_len*4, 11 + 2*V(tf_subgraph)$out_degree)
      
      # ---- Saving Plot ----
      tiff(paste0(name, "_TF_Hierarchy_clean 7.tiff"),
           res = 1200,
           width = 16,
           height = 16,
           units = "in",
           compression = "lzw")
      
      par(mar = c(1, 1, 3, 1))
      set.seed(123)
      plot(tf_subgraph,
           layout = lay,
           rescale = TRUE,
           vertex.label = V(tf_subgraph)$name,
           vertex.label.cex = 1.5,
           vertex.size =V(tf_subgraph)$size,
           vertex.label.color = "black",   # inside circle
           vertex.label.family = "sans",
           vertex.frame.color = "black",
           edge.arrow.size = 1.2,
           edge.width = 1.4,
           edge.color = "gray40",
           edge.curved = 0.2,
           main = paste("TF Regulatory Hierarchy -", name))
      legend("bottomright",
             legend = c("Layer 1 (Upstream)",
                        "Layer 2 (Intermediate)",
                        "Layer 3 (Downstream)"),
             col = layer_colors,
             pch = 19,
             ncol = 3, 
             cex = 0.9,
             bty = "n")
      dev.off()
      
    } else {
      cat("No TF-TF regulatory relationships in", name, "\n")
    }
    
    # ---- Centrality ----
    V(g)$betweenness <- betweenness(g, directed = TRUE)
    V(g)$closeness   <- closeness(g, mode = "all")
    node_centrality <- data.frame(gene        = V(g)$name,
                                  type        = V(g)$type,
                                  degree      = V(g)$degree,
                                  in_degree   = V(g)$in_degree,
                                  out_degree  = V(g)$out_degree,
                                  betweenness = V(g)$betweenness,
                                  closeness   = V(g)$closeness)
    
    top_central <- node_centrality %>%
                   arrange(desc(betweenness)) %>%
                   head(10)
    
    return(list(tf_subgraph = tf_subgraph,
                node_centrality = node_centrality,
                top_central = top_central))
    
  }) 
  names(Subnetworks) <- names(Networks)
  
}



# 9: Validation against known TF-target relationships
#-------------------------------------------------------------------------------
{
  # Using DoRothEA database which contains curated TF-target interactions
  # with confidence levels (A = highest, E = lowest)
  # Filter for high-confidence interactions (levels A, B, and C)
  known_interactions <- dorothea_data %>%
                        filter(confidence %in% c("A", "B", "C")) %>%
                        dplyr::select(tf, target, confidence)
  
  # Overlap of TFs & target between GENIE3 and DoRothEA
  genie3_tfs <- unique(Links$regulatorSymbol)
  genie3_targets <- unique(Links$targetSymbol)
  
  # Filter known interactions to only those with TFs and targets in our dataset
  known_interactions_filtered <- known_interactions %>%
                                 filter(tf %in% genie3_tfs & target %in% genie3_targets)
  known_pairs <- paste(known_interactions_filtered$tf, known_interactions_filtered$target, sep = "_")
  known_interactions_filtered$tf_target_pair <- paste(known_interactions_filtered$tf,
                                                      known_interactions_filtered$target,
                                                      sep = "_")
  
  # Compiling the high confidence links into list
  AllLinks <- Highlinks
  H_links <- Cond.Links$Healthy$high_conf_links
  L_links <- Cond.Links$Lesional$high_conf_links
  NL_links <- Cond.Links$Nonlesional$high_conf_links
  linklist_full <- list(All = AllLinks,
                          H = H_links,
                          NL = NL_links,
                          L = L_links)
  
  Valid_links <- lapply(names(linklist_full), function(x){
      
      df <- linklist_full[[x]]
      cat("Calculating for", x, "\n")
      
      # Generating TF-target pair indentifiers to map to known_interactions
      genie3_pairs <- paste(df$regulatorSymbol,df$targetSymbol, sep = "_")
      
      # GENIE3-predicted interactions supported by the known network
      overlap_pairs <- intersect(genie3_pairs, known_pairs)
      if (length(overlap_pairs) == 0) return(NULL)
      
      # Estimating precision as the proportion of predicted links supported by 
      # known regulatory interactions
      precision <- length(overlap_pairs) / length(genie3_pairs)
      cat("Precision:", round(precision * 100, 2), "%\n\n")

      validated_links <- df[genie3_pairs %in% overlap_pairs, ]
      validated_links$tf_target_pair <- paste(validated_links$regulatorSymbol,validated_links$targetSymbol,sep = "_")
      validated_links <- merge(validated_links,
                               known_interactions_filtered[, c("tf_target_pair", "confidence")],
                               by = "tf_target_pair")
      validated_links %>% arrange(desc(weight))
      
    })
  names(Valid_links) <- names(linklist_full)

  # Combining validated links
  all_validated <- bind_rows(Valid_links)
  all_validated$Condition <- ifelse(all_validated$RegulatoryPair  %in% H_pair$RegulatoryPair, "Healthy", 
                                    ifelse(all_validated$RegulatoryPair  %in% L_pair$RegulatoryPair, "Lesional",
                                           ifelse(all_validated$RegulatoryPair  %in% NL_pair$RegulatoryPair, "Non-Lesional", 
                                                  "Collated")))
  
  # Exporting all validated predictions
  write_xlsx(all_validated, "all validated_tf_target_predictions.xlsx")
  write_xlsx(all_validated %>% filter(!Condition == "Collated"), "condition validated_tf_target_predictions.xlsx")
  
}



# STEP 10: Extracting regulatory Network of 5 Transcription factors shared by 
# Healthy and Non-lesional skin and SOX17
#-------------------------------------------------------------------------------
{
  # Compiling the high confidence links into list
  AllLinks <- Highlinks
  H_links <- Cond.Links$Healthy$high_conf_links
  L_links <- Cond.Links$Lesional$high_conf_links
  NL_links <- Cond.Links$Nonlesional$high_conf_links
  linklist_full <- list(All = AllLinks,
                        H = H_links,
                        NL = NL_links,
                        L = L_links)
  
  # Extracting the regulatory links of the 5 TFs
  # The network is visualised using Cytoscape ~~~
  HNtf_names <- c("POU4F1", "HOXA6", "INSM1","ZNF93", "SOX30")
  HNtf <- lapply(linklist_full, function(df) {
    df <- df %>%
          dplyr::filter( regulatorSymbol %in% HNtf_names | targetSymbol %in% HNtf_names)
    df$targetSymbol <- ifelse(is.na(df$targetSymbol), df$targetGene, df$targetSymbol )
    return(df)
    })
  writexl::write_xlsx(HNtf, "Regulatory network of HN TFs.xlsx")
  # The network is visualised using Cytoscape ~~~
  
  # Extracting the regulatory links of the SOX17
  # The network is visualised using Cytoscape ~~~
  SOX17 <- lapply(linklist_full, function(df) {
    df <- df %>%
          dplyr::filter( regulatorSymbol == "SOX17" | targetSymbol == "SOX17")
    df$targetSymbol <- ifelse(is.na(df$targetSymbol), df$targetGene, df$targetSymbol )
    return(df)
    })
  writexl::write_xlsx(SOX17, paste0("Regulatory network of SOX17.xlsx"))
  
  }






