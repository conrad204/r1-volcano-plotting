# volcanoplot.R
# Standalone R script for volcano plot generation from differential expression data

library(ggplot2)
library(ggrepel)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

# ---- Data Loading ----

# Load differential expression data from a CSV file
# Arguments:
#   file_path: Path to CSV file with columns: EnsemblID, log2FoldChange, padj
# Returns: Data frame with expression data
load_expression_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("CSV file not found: %s", file_path))
  }
  
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Validate required columns

  required_cols <- c("EnsemblID", "log2FoldChange", "padj")
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }
  
  data
}

# ---- Ensembl ID Mapping ----

# Map Ensembl IDs to gene symbols using mouse annotation database
# Arguments:
#   data: Data frame with EnsemblID column
# Returns: Data frame with added SYMBOL column
map_ensembl_to_symbols <- function(data) {
  # Remove version suffix from Ensembl IDs (e.g., ENSMUSG...1.3 -> ENSMUSG...1)
  data$EnsemblID_NoVersion <- gsub("\\..*", "", data$EnsemblID)
  
  # Check if GeneName column exists, otherwise map Ensembl -> SYMBOL

  if ("GeneName" %in% names(data)) {
    data$SYMBOL <- data$GeneName
    message("Using existing GeneName column for gene symbols")
  } else {
    ensembl_ids <- unique(data$EnsemblID_NoVersion)
    gene_symbols_df <- AnnotationDbi::select(
      org.Mm.eg.db,
      keys = ensembl_ids,
      columns = "SYMBOL",
      keytype = "ENSEMBL"
    )
    
    data <- merge(data, gene_symbols_df,
                  by.x = "EnsemblID_NoVersion", by.y = "ENSEMBL", all.x = TRUE)
    message("Mapped Ensembl IDs to gene symbols")
  }
  
  # Precompute -log10(padj)
  data$minusLog10Padj <- -log10(data$padj)
  
  data
}

# Load and prepare expression data (convenience wrapper)
# Arguments:
#   file_path: Path to CSV file
# Returns: Prepared data frame with gene symbols and -log10(padj)
prepare_expression_data <- function(file_path) {
  data <- load_expression_data(file_path)
  data <- map_ensembl_to_symbols(data)
  data
}

# ---- Interesting Genes Management ----

# Create a new interesting genes list
# Arguments:
#   genes: Optional character vector of initial genes
# Returns: Character vector of interesting genes
create_interesting_genes <- function(genes = character(0)) {
  unique(trimws(genes))
}

# Add genes to an interesting genes list
# Arguments:
#   interesting_genes: Current character vector of interesting genes
#   genes: Character vector of genes to add
#   data: Optional data frame to validate genes exist in SYMBOL column
# Returns: Updated character vector of interesting genes
add_interesting_genes <- function(interesting_genes, genes, data = NULL) {
  genes <- trimws(genes)
  
  if (!is.null(data)) {
    valid_genes <- genes[genes %in% data$SYMBOL]
    invalid_genes <- setdiff(genes, valid_genes)
    if (length(invalid_genes) > 0) {
      message(sprintf("Genes not found in data: %s", paste(invalid_genes, collapse = ", ")))
    }
    genes <- valid_genes
  }
  
  unique(c(interesting_genes, genes))
}

# Remove genes from an interesting genes list
# Arguments:
#   interesting_genes: Current character vector of interesting genes
#   genes: Character vector of genes to remove
# Returns: Updated character vector of interesting genes
remove_interesting_genes <- function(interesting_genes, genes) {
  setdiff(interesting_genes, trimws(genes))
}

# Check if a gene exists in the data
# Arguments:
#   gene: Gene symbol to check
#   data: Data frame with SYMBOL column
# Returns: Logical indicating if gene exists
gene_exists <- function(gene, data) {

  trimws(gene) %in% data$SYMBOL
}

# ---- Summary Statistics ----

# Get summary counts for expression data
# Arguments:
#   data: Data frame with expression data
#   lfc: Log2 fold change threshold (default 1.5)
#   p_cut: Adjusted p-value cutoff (default 0.05)
#   interesting_genes: Character vector of interesting genes
# Returns: List with summary counts
get_summary <- function(data, lfc = 1.5, p_cut = 0.05, interesting_genes = character(0)) {
  list(
    n_total = nrow(data),
    n_sig_up = sum(data$log2FoldChange > lfc & data$padj < p_cut, na.rm = TRUE),
    n_sig_down = sum(data$log2FoldChange < -lfc & data$padj < p_cut, na.rm = TRUE),
    n_interesting = sum(data$SYMBOL %in% interesting_genes, na.rm = TRUE)
  )
}

# ---- Volcano Plot Generation ----

# Generate a volcano plot
# Arguments:
#   data: Data frame with columns: log2FoldChange, padj (or minusLog10Padj), SYMBOL
#   interesting_genes: Character vector of genes to highlight in green
#   lfc: Log2 fold change threshold for significance (default 1.5)
#   minlogp: -log10(padj) threshold for significance (default -log10(0.05))
#   show_threshold_labels: Logical; show labels for genes passing thresholds (default TRUE)
#   title: Plot title (default "Volcano Plot")
# Returns: ggplot object
generate_volcano_plot <- function(data,
                                   interesting_genes = character(0),
                                   lfc = 1.5,
                                   minlogp = -log10(0.05),
                                   show_threshold_labels = TRUE,
                                   title = "Volcano Plot") {
  
 df <- data
  
  # Ensure minusLog10Padj exists
  if (is.null(df$minusLog10Padj)) {
    df$minusLog10Padj <- -log10(df$padj)
  }
  
  # Assign dot colors based on thresholds
  # Red: log2FC > lfc AND -log10(padj) > minlogp
  # Blue: log2FC < -lfc AND -log10(padj) > minlogp
  # Green: interesting genes
  df$dot_color <- "gray"
  df$dot_color[!is.na(df$log2FoldChange) & df$log2FoldChange > lfc & 
               !is.na(df$minusLog10Padj) & df$minusLog10Padj > minlogp] <- "red"
  df$dot_color[!is.na(df$log2FoldChange) & df$log2FoldChange < -lfc & 
               !is.na(df$minusLog10Padj) & df$minusLog10Padj > minlogp] <- "blue"
  df$dot_color[df$SYMBOL %in% interesting_genes] <- "green"
  
  # Determine which points get labels
  threshold_labels <- (!is.na(df$log2FoldChange) & (abs(df$log2FoldChange) >= lfc)) &
    (!is.na(df$minusLog10Padj) & (df$minusLog10Padj >= minlogp))
  
  df$label <- ifelse(
    ((threshold_labels & show_threshold_labels) | (df$SYMBOL %in% interesting_genes)),
    ifelse(is.na(df$SYMBOL), df$EnsemblID, df$SYMBOL),
    NA
  )
  
  # Calculate dynamic axis limits
  max_abs_lfc <- max(abs(df$log2FoldChange), na.rm = TRUE)
  x_limit <- max_abs_lfc + 0.5
  
  max_y <- max(df$minusLog10Padj, na.rm = TRUE)
  y_limit <- max_y + 3
  
  # Build plot with layered points (gray first, then colored, green on top)
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
    labs(title = title,
         x = "log2(Fold Change)",
         y = "-log10(Adjusted p-value)")
  
  p
}

# Save volcano plot to file
# Arguments:
#   plot: ggplot object from generate_volcano_plot
#   file_path: Output file path (supports png, pdf, etc.)
#   width: Plot width in inches (default 10)
#   height: Plot height in inches (default 8)
#   dpi: Resolution for raster formats (default 300)
save_volcano_plot <- function(plot, file_path, width = 10, height = 8, dpi = 300) {
  ggsave(file_path, plot = plot, width = width, height = height, dpi = dpi)
  message(sprintf("Plot saved to: %s", file_path))
}
