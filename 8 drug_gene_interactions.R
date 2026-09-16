#======================================================================================
#                 DRUG–GENE INTERACTION API PIPELINE
# DrugCentral, ChEMBL, OpenTargets, PharmaGKB, KEGG DBGET, DrugBank, HCDT, TTD, DGidb
#======================================================================================
 
# Drug–gene interaction information for a predefined set of target genes from
# DrugCentral and ChEMBL. 
# Interaction records containing the target gene, interacting drug, 
# interaction type, evidence/source information, publication identifiers, 
# database URLs, and activity values. 



# Packages
{
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(readxl)
  library(xml2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(writexl)
  library(data.table)
  library(ggplot2)
  library(scales)
  
}



# Utitlity functions
#-------------------------------------------------------------------------------
{
  # Operator
  # Defines a null-coalescing operator to return a fallback value when the primary 
  # value is NULL or empty.
  `%||%` <- function(a,b) if(!is.null(a) && length(a)>0) a else b
  
  
  
  # JSON
  # Parses JSON responses returned by web APIs. 
  # The function supports common R HTTP response object formats and returns NULL 
  # when a valid JSON response cannot be extracted.
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
  
  
  
  # Executes single-gene database queries 
  # Database-specific query functions while preventing an error by one gene from  
  # interrupting the complete drug–gene interaction analysis.
  safe_query <- function(fn, gene) {
    tryCatch(
      fn(gene),
      error = function(e) {
        message("[ERROR] ", gene, " in ", deparse(substitute(fn)), ": ", conditionMessage(e))
        tibble()   
      }
    )
  }
  
  
  
  # Organizes database-specific query results 
  # Applies a database-specific query function to each gene and store the resulting 
  # interaction records as a named list, with one element corresponding to each gene.
  build_db_list <- function(db_name, fn, genes) {
    results <- setNames(
      lapply(genes, function(g) safe_query(fn, g)),
      genes
    )
    return(results)
  }
  
  
  
  # Drug–gene interaction extraction pipeline 
  # Executes all database-specific query functions for the predefined genes of 
  # interest and organizes the results into a nested database-by-gene list.
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
  
}



# Functions to query Databases
#-------------------------------------------------------------------------------
{
  # 1) Drug Central
  #-----------------------------------------------------------------------------
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

  
  # 2) ChEMBL 
  #-----------------------------------------------------------------------------
  # Query ChEMBL for drug–target activity records associated with each gene
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
  
  # Applies the ChEMBL extraction function to the complete set of target genes.
  query_chembl <- function(genes) {
    bind_rows(lapply(genes, query_chembl_single))
  }

  
  # 3) OpenTargets
  #-----------------------------------------------------------------------------
  query_opentargets <- function(gene_symbol, size = Inf) {
    
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

  
  # 4) PharmaGKB
  #-----------------------------------------------------------------------------
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

  
  # 5) KEGG DBGET
  #-----------------------------------------------------------------------------
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

  
  # 6) DrugBank 
  #-----------------------------------------------------------------------------
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

  
  # 7) HCDT
  #-----------------------------------------------------------------------------
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
  
  
  # 8) TTD
  #-----------------------------------------------------------------------------
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

}



# Functions to perform Drug Gene Interaction analysis
#-------------------------------------------------------------------------------
{
  # Harmonization of drug–gene interaction results 
  # Integrates drug–gene interaction records from multiple databases into a 
  # standardized master dataset containing common gene, drug, source, regulatory, 
  # interaction, and interaction-score fields.
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
  
  
  # Drug evidence summary 
  # Identifies approved drugs and calculates a standardized binding score for 
  # drug–gene interactions with available chemical interaction values.
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
    
  }
  
  
  
  # Binding score distribution and threshold analysis 
  # Evaluate the distribution of chemical interaction scores, identify 
  # percentile-based score thresholds, and quantify the number of interactions 
  # exceeding each threshold.
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

}



# Extraction and Analysis
#-------------------------------------------------------------------------------
{
  # Defining genes to study interactions of
  tg <- c("DLGAP5",	"NCAPG",	"CHI3L2",	"SPRR2E",	"LYPD2", "OAS2", "PARP9", "OAS1", "IFI44", "MX1",	
          "CCL13",	"GALNT6",	"POLR3G",	"TTC39A",	"CCNB1",	"NELL2", "PDZK1",	"INSIG1",	"FAR2",	
          "AADACL3","CTNNA1",	"PGK1",	"UBE2N",	"PGAM1", "SERPINA12",	"CLDN1",	"FLG2",	"SCEL")
  
  # Retrieving drug–gene interaction records for the target hub genes
  # from the integrated drug interaction resources.
  TgDrugs <- dgInt(tg)
  
  # Loading literature-supported drug–gene interaction records from DGIdb
  # for complementary evidence and interaction-score annotation and Incorporating
  # DGIdb evidence.
  TargetDGIdbres <-  read_delim("Targets gene_interaction_results-15_06_2026.tsv")
  TargetDGIdbres$Source <- rep("DGidb", nrow(TargetDGIdbres))
  range(TargetDGIdbres$`interaction score`)
  TargetDGIdbres <- TargetDGIdbres %>%
    dplyr::rename(DGIdb_Interaction_Score = `interaction score`)
  
  
  # Integrating and harmonizing database-derived interactions along with DGIdb 
  # evidence into a standardized master drug–gene interaction dataset.
  TgDrugs <- harmonize_drug_results( Drugs_out = TgDrugs,
                                     DGIdbres = TargetDGIdbres,
                                     save_file = "Target_Drugs.xlsx")
  
  # Identifying the approved drugs and calculate standardized chemical interaction
  # scores for the target gene–drug interactions.
  DrugEvidence_tg <- drug_evidence_summary2(Drugs = TgDrugs$Master_Drugs,
                                            save_prefix = "Target_EV")
  
  
  # Analyzing chemical interaction scores distribution and determining
  # percentile-based thresholds for prioritizing stronger interactions.
  BindingAnalysis_tg <- analyze_binding_scores(DrugEvidence_tg$BindingScore, 
                                               score_col = "Log10ChEMBLintscore",
                                               probs = c(0.25, 0.50, 0.75, 0.95, 0.97, 0.99),
                                               make_plot = TRUE)
  # $density
  # Call:
  #  density.default(x = na.omit(intscore))
  # Data: na.omit(intscore) (10 obs.);	Bandwidth 'bw' = 1.427
  # x                  y            
  # Min.   :-0.09114   Min.   :0.0004345  
  # 1st Qu.: 3.60141   1st Qu.:0.0154052  
  # Median : 7.29397   Median :0.0850165  
  # Mean   : 7.29397   Mean   :0.0675350  
  # 3rd Qu.:10.98653   3rd Qu.:0.1068076  
  # Max.   :14.67908   Max.   :0.1395308 
  
  # Quantiles and no. of interactions at each quantiles
  # 25%       50%       75%       95%       97%       99% 
  # 4.602932  8.037860  9.136797  10.139626 10.242952 10.346277
  # 7         5         3         1         1         1 
  
  # Inspecting genes and drugs retained after score filtering
  unique( BindingAnalysis_tg$Cleaned_Drugs$Gene )
  # [1] "PARP9" 
  unique(BindingAnalysis_tg$Cleaned_Drugs$Drug)
  # [1] "CHEMBL4287262" "CHEMBL449635"  "CHEMBL4287655" "CHEMBL4208737" "CHEMBL1438938" "ATAMPARIB"     "RUCAPARIB"     "CHEMBL1488758"
  
  
  # Retaining one record per gene–drug pair to obtain a non-redundant
  # interaction dataset for downstream interpretation.
  BindingAnalysis_tg$Cleaned_Drugs <- BindingAnalysis_tg$Cleaned_Drugs %>%
    dplyr::rename(`-Log10ChEMBLScore` = Log10ChEMBLintscore)%>%
    dplyr::group_by(Gene)%>%
    dplyr::distinct(Gene, Drug, .keep_all = TRUE)
  write_xlsx( BindingAnalysis_tg$Cleaned_Drugs, "Drug-Gene Interaction Target revised.xlsx")
  
}
