library(Seurat)
library(ggplot2)
library(ggsignif)

is_significant <- function(obj, gene, conditions=NULL, asterisk = F) {
  
  df <- NULL
  if (is.null(conditions)) {
    conditions <- as.vector(unique(Idents(obj)))
  }
  for (con in conditions) {
    cond <- obj[gene, Idents(obj) == con]
    cond_df <- data.frame(
      gene = as.vector(GetAssayData(cond, slot = 'data')), 
      group = con, 
      split = con)
    rownames(cond_df) <- colnames(cond)
    df <- rbind(df, cond_df)
  }
  
  p <- ggplot(data=df, aes(x=factor(group, level = conditions), y=gene)) + 
    geom_violin(aes(fill=group)) +
    geom_jitter(height = 0, width = 0.1) +
    geom_signif(comparisons = combn(conditions, 2, simplify = F),
                map_signif_level = asterisk, textsize=5, test = 'wilcox.test',
                step_increase = .1) +
    labs(x = gene, y = '')
  
  return (p)
}

get_expr_data <- function(seu, gene, data_name) {
  df <- NULL
  seu_split <- SplitObject(seu[gene, ])
  for (seu in seu_split) {
    ident <- unique(Idents(seu))
    seu_df <- data.frame(expr = as.vector(Seurat::GetAssayData(seu, slot = "data")), group=ident, gene=gene, data=data_name)
    df <- rbind(df, seu_df)
  }
  
  return (df)
}

plot_violins <- function(fung, li, genes, show_points=T, sig=F) {
  for (gene in genes) {
    li_df <- get_expr_data(li, gene, 'Li et al., 2020')
    fung_df <- get_expr_data(fung, gene, 'In vitro')
    df <- rbind(fung_df, li_df)
    df$group <- factor(df$group, levels = c("hADE", "VFG", "hFG", "hAL", "hMG", "hHG"))
    df$data <- factor(df$data, levels = c("In vitro", "Li et al., 2020"))
    p <- ggplot(df, aes(x=group, y=expr, fill=group)) + 
      geom_violin(scale = "width") + 
      scale_fill_manual(values=c("#797979", "#4E8F32", "#025393", "#25662C", "#612E7F", "#931813")) +
      labs(x = gene, y = 'Expression', fill = 'Cell types') +
      theme_classic()
    
    if (show_points) {
      p <- p + geom_jitter(height = 0, width = 0.1)
    }
    
    if (!sig) {
      p <- p + facet_grid(~data, scales = "free", space = "free")
    } else {
      combinations <- c(
        utils::combn(as.vector(unique(Idents(li))), 2, simplify = F),
        utils::combn(as.vector(unique(Idents(fung))), 2, simplify = F)
      )
      p <- p + geom_signif(comparisons = combinations, test = "wilcox.test", 
                           step_increase = 0.1, size = 0.5, textsize = 3, 
                           vjust = 0, map_signif_level = F)
    }
    
    plot(p)
  }
}