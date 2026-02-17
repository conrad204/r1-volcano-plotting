# plumber.R (top of file)

options(
  plumber.serializer.png.width  = 4000,
  plumber.serializer.png.height = 3000,
  plumber.serializer.png.res    = 200
)

library(plumber)
library(jsonlite)
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Try several candidate CSV locations (RES_FILE env, container /data, repo-local data/)
candidates <- c(
  Sys.getenv("RES_FILE", ""),
  "/data/example.csv",
  file.path(getwd(), "data", "example.csv"),
  "data/example.csv"
)
candidates <- unique(candidates[nzchar(candidates)])
existing <- candidates[file.exists(candidates)]

if (length(existing) == 0) {
  stop(sprintf("CSV not found. Looked in: %s\nSet RES_FILE or place example.csv in the project directory.",
               paste(candidates, collapse = ", ")))
}

res_file <- existing[1]
message(sprintf("Using CSV: %s", res_file))
res <- read.csv(res_file, stringsAsFactors = FALSE)



# ---- Load & prepare data once ----
# data already loaded above into `res` from the chosen candidate path

# Remove version suffix from Ensembl IDs (e.g., ENSMUSG...1.3 -> ENSMUSG...1)
res$EnsemblID_NoVersion <- gsub("\\..*", "", res$EnsemblID)

# Check if GeneName column exists, otherwise map Ensembl -> SYMBOL
if ("GeneName" %in% names(res)) {
  # Use existing GeneName column as SYMBOL
  res$SYMBOL <- res$GeneName
  message("Using existing GeneName column for gene symbols")
} else {
  # Map Ensembl -> SYMBOL (mouse)
  ensembl_ids <- unique(res$EnsemblID_NoVersion)
  gene_symbols_df <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = ensembl_ids,
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
  
  # Left join to add SYMBOL
  res <- merge(res, gene_symbols_df,
               by.x = "EnsemblID_NoVersion", by.y = "ENSEMBL", all.x = TRUE)
  message("Mapped Ensembl IDs to gene symbols")
}

# Precompute -log10(padj)
res$minusLog10Padj <- -log10(res$padj)

# ---- State (lives across requests) ----
state <- new.env(parent = emptyenv())
state$interesting_genes <- c()

# Utility: normalize a gene symbol (trim + keep original case comparison exact)
.crumb <- function(x) trimws(x)

# Recompute colors/labels whenever interesting_genes changes
recompute_annotations <- function() {
  df <- res

  # dot colors
  df$dot_color <- "gray"
  df$dot_color[df$log2FoldChange > 1.5 & df$padj < 0.05] <- "red"
  df$dot_color[df$log2FoldChange < -1.5 & df$padj < 0.05] <- "blue"
  df$dot_color[df$SYMBOL %in% state$interesting_genes] <- "green"

  # labels
  df$label <- ifelse(
    (df$minusLog10Padj > 60 |
       (df$minusLog10Padj > 20 & df$log2FoldChange < -4) |
       df$SYMBOL %in% state$interesting_genes),
    ifelse(is.na(df$SYMBOL), df$EnsemblID, df$SYMBOL),
    NA
  )

  state$res_annotated <- df
  invisible(TRUE)
}
recompute_annotations()

# ---- CORS so you can call from a file:// HTML page or localhost frontend ----
#* @filter cors
function(req, res){
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, DELETE, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type")
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$status <- 200
    return(list())
  }
  plumber::forward()
}

# ---- Helpers ----
gene_exists_in_res <- function(gene) {
  # exact match on SYMBOL
  gene %in% res$SYMBOL
}

# ---- API: list interesting genes ----
#* Get current interesting genes
#* @get /genes
function(){
  list(interesting_genes = state$interesting_genes,
       count = length(state$interesting_genes))
}

# ---- API: add a gene if it exists in res$SYMBOL ----
# Accepts query (?gene=Notch3) or form-urlencoded or JSON { "gene": "Notch3" }

#* Add gene to interesting list (if present in res$SYMBOL)
#* @post /genes
function(req, gene = NULL){
  # attempt to read JSON body if no form/query param
  if (is.null(gene)) {
    if (!is.null(req$postBody) && nzchar(req$postBody)) {
      body <- tryCatch(jsonlite::fromJSON(req$postBody), error = function(e) NULL)
      if (!is.null(body$gene)) gene <- body$gene
    }
  }
  if (is.null(gene) || !nzchar(gene)) {
    return(list(ok = FALSE, message = "Provide 'gene'"))
  }

  gene <- .crumb(gene)

  if (!gene_exists_in_res(gene)) {
    return(list(ok = TRUE, found = FALSE,
                message = sprintf("'%s' not found in res$SYMBOL", gene)))
  }

  if (gene %in% state$interesting_genes) {
    return(list(ok = TRUE, found = TRUE, already_present = TRUE,
                message = sprintf("'%s' already in interesting_genes", gene),
                interesting_genes = state$interesting_genes))
  }

  state$interesting_genes <- c(state$interesting_genes, gene)
  state$interesting_genes <- unique(state$interesting_genes)
  recompute_annotations()

  list(ok = TRUE, found = TRUE, added = gene,
       interesting_genes = state$interesting_genes)
}

# ---- API: remove a gene from interesting list ----
#* Remove gene from interesting list
#* @delete /genes
function(gene){
  gene <- .crumb(gene)
  if (!gene %in% state$interesting_genes) {
    return(list(ok = TRUE, removed = FALSE,
                message = sprintf("'%s' not in interesting_genes", gene),
                interesting_genes = state$interesting_genes))
  }
  state$interesting_genes <- setdiff(state$interesting_genes, gene)
  recompute_annotations()
  list(ok = TRUE, removed = gene, interesting_genes = state$interesting_genes)
}

# ---- API: reset/clear the interesting genes ----
#* Reset (clear) the interesting_genes list
#* @post /genes/reset
function(){
  state$interesting_genes <- character(0)
  recompute_annotations()
  list(ok = TRUE, cleared = TRUE, interesting_genes = state$interesting_genes)
}

# ---- API: upload custom CSV ----
#* Upload a custom CSV file to replace the current dataset
#* @post /upload
#* @parser csv
function(req){
  tryCatch({
    # The CSV parser automatically reads the uploaded CSV into a data frame
    uploaded_data <- req$body
    
    # validate we got data
    if (is.null(uploaded_data) || nrow(uploaded_data) == 0) {
      return(list(ok = FALSE, message = "No data in uploaded CSV"))
    }
    
    # basic validation - check for required columns
    required_cols <- c("EnsemblID", "log2FoldChange", "padj")
    missing <- setdiff(required_cols, names(uploaded_data))
    if (length(missing) > 0) {
      return(list(ok = FALSE, 
                  message = sprintf("Missing required columns: %s", paste(missing, collapse = ", "))))
    }
    
    # process the uploaded data (same as initialization)
    uploaded_data$EnsemblID_NoVersion <- gsub("\\..*", "", uploaded_data$EnsemblID)
    
    # Check if GeneName column exists, otherwise map Ensembl -> SYMBOL
    if ("GeneName" %in% names(uploaded_data)) {
      # Use existing GeneName column as SYMBOL
      uploaded_data$SYMBOL <- uploaded_data$GeneName
      message("Using existing GeneName column for gene symbols")
    } else {
      # Map Ensembl -> SYMBOL
      ensembl_ids_upload <- unique(uploaded_data$EnsemblID_NoVersion)
      gene_symbols_upload <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = ensembl_ids_upload,
        columns = "SYMBOL",
        keytype = "ENSEMBL"
      )
      
      # merge symbols
      uploaded_data <- merge(uploaded_data, gene_symbols_upload,
                            by.x = "EnsemblID_NoVersion", by.y = "ENSEMBL", all.x = TRUE)
      message("Mapped Ensembl IDs to gene symbols")
    }
    
    # compute -log10(padj)
    uploaded_data$minusLog10Padj <- -log10(uploaded_data$padj)
    
    # replace global res with uploaded data
    res <<- uploaded_data
    
    # clear interesting genes list (new dataset)
    state$interesting_genes <- character(0)
    
    # recompute annotations with new data
    recompute_annotations()
    
    list(ok = TRUE, 
         message = "CSV uploaded and processed successfully",
         rows = nrow(res),
         columns = ncol(res))
         
  }, error = function(e) {
    list(ok = FALSE, message = sprintf("Error processing CSV: %s", e$message))
  })
}

# ---- API: check membership only ----
#* Check if a gene exists in res$SYMBOL and/or is interesting
#* @get /check
function(gene){
  gene <- .crumb(gene)
  list(
    query = gene,
    exists_in_res = gene_exists_in_res(gene),
    in_interesting = gene %in% state$interesting_genes
  )
}

# ---- API: return a quick summary (counts) ----
#* Summary counts for current annotation
#* @get /summary
function(){
  df <- state$res_annotated
  list(
    n_total = nrow(df),
    n_sig_up = sum(df$log2FoldChange > 1.5 & df$padj < 0.05, na.rm = TRUE),
    n_sig_down = sum(df$log2FoldChange < -1.5 & df$padj < 0.05, na.rm = TRUE),
    n_interesting = sum(df$SYMBOL %in% state$interesting_genes, na.rm = TRUE)
  )
}

# ---- API: serve volcano plot as PNG ----
#* @serializer png
#* @get /plot
function(lfc = 1.5, minlogp = -log10(0.05), p_cut = 0.05, show_threshold_labels = TRUE){
  # accept thresholds from query params (lfc = log2 fold-change threshold,
  # minlogp = -log10(padj) threshold). p_cut is an underlying padj cutoff used to color points.
  lfc <- as.numeric(lfc)
  minlogp <- as.numeric(minlogp)
  p_cut <- as.numeric(p_cut)

  # interpret show_threshold_labels - accept strings 'false'/'0' as FALSE
  if (is.character(show_threshold_labels)){
    stl <- tolower(show_threshold_labels)
    if (stl %in% c('false','0','f','no','n')) show_threshold_labels <- FALSE
    else show_threshold_labels <- TRUE
  } else {
    show_threshold_labels <- as.logical(show_threshold_labels)
  }

  df <- res
  if (is.null(df$minusLog10Padj)) df$minusLog10Padj <- -log10(df$padj)

  # dot colors based on thresholds
  # Red: log2FC > lfc AND -log10(padj) > minlogp
  # Blue: log2FC < -lfc AND -log10(padj) > minlogp
  df$dot_color <- "gray"
  df$dot_color[!is.na(df$log2FoldChange) & df$log2FoldChange > lfc & 
               !is.na(df$minusLog10Padj) & df$minusLog10Padj > minlogp] <- "red"
  df$dot_color[!is.na(df$log2FoldChange) & df$log2FoldChange < -lfc & 
               !is.na(df$minusLog10Padj) & df$minusLog10Padj > minlogp] <- "blue"
  df$dot_color[df$SYMBOL %in% state$interesting_genes] <- "green"

  # threshold-based labels (boolean vector)
  threshold_labels <- (!is.na(df$log2FoldChange) & (abs(df$log2FoldChange) >= lfc)) &
    (!is.na(df$minusLog10Padj) & (df$minusLog10Padj >= minlogp))

  # final label if (threshold label and allowed) OR in interesting genes
  df$label <- ifelse(
    ((threshold_labels & show_threshold_labels) | (df$SYMBOL %in% state$interesting_genes)),
    ifelse(is.na(df$SYMBOL), df$EnsemblID, df$SYMBOL),
    NA
  )

  # Calculate dynamic axis limits
  # For x-axis: find max absolute value, then add 0.5 and make symmetric
  max_abs_lfc <- max(abs(df$log2FoldChange), na.rm = TRUE)
  x_limit <- max_abs_lfc + 0.5
  
  # For y-axis: find max value and add 3
  max_y <- max(df$minusLog10Padj, na.rm = TRUE)
  y_limit <- max_y + 3
  
  p <- ggplot(df, aes(x = log2FoldChange, y = minusLog10Padj)) +
    geom_point(data = subset(df, dot_color == "gray"),
               aes(color = dot_color), alpha = 0.7) +
    geom_point(data = subset(df, dot_color == "blue"),
               aes(color = dot_color), alpha = 0.7) +
    geom_point(data = subset(df, dot_color == "red"),
               aes(color = dot_color), alpha = 0.7) +
    geom_point(data = subset(df, dot_color == "green"),
               aes(color = dot_color), alpha = 1) +
    scale_color_identity() +
    geom_text_repel(
      data = subset(df, !is.na(label)),
      aes(label = label),
      color = "black",
      size = 3,
      fontface = "bold.italic",
      segment.color = "black",
      segment.size = 0.4,
      min.segment.length = 0.1,
      max.overlaps = Inf
    ) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dotted", color = "gray50", linewidth = 0.8) +
    geom_hline(yintercept = minlogp, linetype = "dotted", color = "gray50", linewidth = 0.8) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    scale_y_continuous(limits = c(0, y_limit)) +
    theme_bw() +
    labs(title = "Volcano Plot",
         x = "log2(Fold Change)",
         y = "-log10(Adjusted p-value)")
  print(p)
}
