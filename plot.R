##R4.3.1
#' Generate a Bidirectional EWAS Manhattan Plot
#'
#' @param df Data frame containing EWAS results (columns: Chromosome, Position, p, beta, cpgs)
#' @param thresholds Numeric vector (length 2) defining significance thresholds (default: c(9e-8, 1e-5))
#' @param colors Named list for plot colors (e.g., list(threshold1="red", threshold2="blue", points="grey60", highlight="red"))
#' @param highlight_file CSV file with SNPs to highlight (column: `cpgs`)
#' @param config_file YAML file for custom settings
#' @return A ggplot2 object
#' @export
BD_manhattan <- function(df, thresholds = c(9e-8, 1e-5), colors = NULL, highlight_file = NULL, config_file = NULL) {
  
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(yaml)
  
  if (!is.null(config_file)) {
    config <- yaml::read_yaml(config_file)
    if (!is.null(config$thresholds)) thresholds <- config$thresholds
    if (!is.null(config$colors)) colors <- config$colors
    if (!is.null(config$highlight_file)) highlight_file <- config$highlight_file
  }
  
  if (is.null(colors)) {
    colors <- list(threshold1 = "red", threshold2 = "blue", points = "grey60", highlight = "red")
  }
  
  if (!is.null(highlight_file)) {
    highlight_snps <- read.csv(highlight_file)$cpgs
  } else {
    highlight_snps <- df %>% filter(p < thresholds[2]) %>% pull(cpgs)
  }
  
  df <- df %>% 
    group_by(Chromosome) %>% 
    mutate(tot = cumsum(max(Position)) - max(Position)) %>%
    ungroup() %>%
    mutate(PositionCum = Position + tot, y = -log10(p) * ifelse(beta > 0, 1, -1))
  
  axisdf <- df %>% group_by(Chromosome) %>% summarize(center = (max(PositionCum) + min(PositionCum)) / 2)
  
  ggplot(df, aes(x = PositionCum, y = y)) +
    geom_point(aes(color = as.factor(Chromosome)), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c(colors$points, "#4197d8"), length(unique(df$Chromosome)) / 2 + 1)) +
    geom_hline(yintercept = -log10(thresholds[2]), linetype = "dashed", color = colors$threshold2) +
    geom_hline(yintercept = log10(thresholds[2]), linetype = "dashed", color = colors$threshold2) +
    geom_hline(yintercept = -log10(thresholds[1]), linetype = "dashed", color = colors$threshold1) +
    geom_hline(yintercept = log10(thresholds[1]), linetype = "dashed", color = colors$threshold1) +
    geom_point(data = subset(df, cpgs %in% highlight_snps), color = colors$highlight, size = 2) +
    geom_text_repel(data = subset(df, p < 1e-5), aes(label = cpgs), size = 4) +
    scale_x_continuous(name = "Chromosome", breaks = axisdf$center, labels = axisdf$Chromosome) +
    scale_y_continuous(name = "-log10(p)", limits = c(-8, 8)) +
    theme_bw() +
    theme(legend.position = "none")
}
