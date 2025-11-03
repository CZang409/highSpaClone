# Download data here: https://drive.google.com/drive/folders/17A0wfkL3YjWkfUko91vZZElBXoz6JW2W?usp=drive_link

##Figure 2
####Fig2b CNV heatmap: see our tutorial website: https://czang409.github.io/highSpaClone/articles/xenium_breast1_code.html
####Fig2d was created by adobe illustrator

##----Figure 2a Method Comparison----
library(ggplot2)
library(dplyr)
library(patchwork)

load('./F2a.RData')

make_cluster_plot <- function(data, colname, title, colors=c("#ebe5c2", "#D57358","#8a508f", "#023047", "#F7CA71")) {
  ggplot(data = data, aes(x = x, y = y, color = factor(.data[[colname]]))) +
    geom_point(size = 0.2) +
    scale_x_reverse() +
    scale_color_manual(name = "cluster", values = colors) +
    theme(
      panel.background  = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_blank(),
      axis.ticks        = element_blank(),
      axis.text         = element_blank(),
      axis.title        = element_blank(),
      plot.title        = element_text(hjust = 0.5, face = 'bold', size = 16),
      plot.title.position = "plot",
      plot.margin       = margin(t = 5, r = 5, b = 5, l = 5),
      panel.border      = element_rect(color = "black", fill = NA)
    ) +
    ggtitle(title) +
    guides(color = "none")
}

p1 <- make_cluster_plot(F2A.data, "highSpaClone", "highSpaClone")
p2 <- make_cluster_plot(F2A.data, "iris", "IRIS")
p3 <- make_cluster_plot(F2A.data, "infercnv", "InferCNV")
p4 <- make_cluster_plot(F2A.data, "truth", "Ground Truth")

# combine four figs
combined_plot <- p4 + p1 + p2 + p3 + plot_layout(ncol = 4)

# save figure
ggsave("./F2 data/Fig2a.png", combined_plot, width = 20, height = 4, dpi = 300)

pdf("./Fig2a.pdf", 20, 4)	 
print (combined_plot)
dev.off()

##----Figure 2c CNV burden----
## total 238 bins
load('./F2c.RData')
custom_colors <- c("CNV_gain" = "#ffadad", "CNV_loss" = "#a0c4ff")
p <- ggplot(Fig2C.data, aes(x = highSpaClone, y = CNV_burdens, fill = CNV_status)) +
  geom_violin(position = position_dodge(width = 0.8),
              color = NA,
              alpha = 0.5,
              width = 0.8,
              trim = TRUE,
              scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6) +  # Remove outliers, reduce box width
  geom_point(data = Fig2C.data %>%
               group_by(highSpaClone, CNV_status) %>%
               sample_n(300),
             aes(x = highSpaClone, y = CNV_burdens, color = CNV_status),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             shape = 16,
             size = 1,
             alpha = 0.5,
             show.legend = FALSE)+
  
  scale_fill_manual(values = custom_colors) +  # Custom colors
  scale_x_discrete(labels = c("DCIS #1", "DCIS #2", "Invasive")) +
  labs(title = "CNV Burdens Across Tumor Clones",
       x = "",
       y = "CNV Burdens",
       fill = "State") +  # Set labels
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # 1. Keep x-axis labels horizontal
    panel.background = element_blank(),  # Remove panel background
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 2. Add black border
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Keep axis lines
    axis.ticks = element_line(color = "black")  # Keep axis ticks
    #plot.background = element_rect(color = "black", size = 1)  # 3. Add black border around the plot
  )+
  coord_cartesian(ylim = c(5, 55))

ggsave("./Fig2c.png", plot = p,
       width = 6.5, height = 4.5, dpi = 300, units = "in")

pdf("./Fig2c.pdf", 6.5, 4.5)	 
print (p)
dev.off()

##----Figure 2e Pseudotime----
load('./Fig2e.RData')
plot_trajectory = function(pseudotime, location, clusterlabels, gridnum, color_in,
                           pointsize = 1, arrowlength = 0.2, arrowsize = 1, 
                           arrowdist = 1, textsize = 22) {
  library(ggplot2)
  library(grid)
  library(dplyr)
  
  pseudotime_use = pseudotime
  info = as.data.frame(location)
  colnames(info) = c("sdimx", "sdimy")
  grids = gridnum
  
  min_x = min(info$sdimx)
  min_y = min(info$sdimy)
  max_x = max(info$sdimx)
  max_y = max(info$sdimy)
  
  x_anchor = seq(min_x, max_x, length.out = grids + 1)
  y_anchor = seq(min_y, max_y, length.out = grids + 1)
  
  count = 0
  start_x_dat = c()
  start_y_dat = c()
  end_x_dat = c()
  end_y_dat = c()
  
  for (num_x in 1:grids) {
    for (num_y in 1:grids) {
      filter_x = which(info$sdimx >= x_anchor[num_x] & info$sdimx <= x_anchor[num_x + 1])
      filter_y = which(info$sdimy >= y_anchor[num_y] & info$sdimy <= y_anchor[num_y + 1])
      points_in_grid = intersect(filter_x, filter_y)
      
      if (length(points_in_grid) > 1 && sum(which.min(pseudotime_use[points_in_grid])) > 0) {
        count = count + 1
        min_point = info[points_in_grid[which.min(pseudotime_use[points_in_grid])], ]
        max_point = info[points_in_grid[which.max(pseudotime_use[points_in_grid])], ]
        start_x_dat[count] = min_point$sdimx
        start_y_dat[count] = min_point$sdimy
        end_x_dat[count] = max_point$sdimx
        end_y_dat[count] = max_point$sdimy
      }
    }
  }
  
  loc1 = info$sdimx
  loc2 = info$sdimy
  time = pseudotime_use + 0.01
  datt = data.frame(time, loc1, loc2)
  datt2 = data.frame(start_x_dat, start_y_dat, end_x_dat, end_y_dat)
  
  datt2$dist <- sqrt((datt2$start_x_dat - datt2$end_x_dat)^2 + (datt2$start_y_dat - datt2$end_y_dat)^2)
  datt2 <- datt2 %>% filter(dist > arrowdist) 
  
  df_plot = info %>%
    mutate(pseudotime = pseudotime_use, cluster = clusterlabels)
  
  pseudotime_arrow_plot <- ggplot() +
    geom_point(data = df_plot %>% filter(is.na(pseudotime)),
               aes(x = sdimx, y = sdimy),
               color = "#eaeaea",
               size = pointsize) +
    geom_point(data = df_plot %>% filter(!is.na(pseudotime)),
               aes(x = sdimx, y = sdimy, color = pseudotime),
               size = pointsize) +
    geom_segment(data = datt2,
                 aes(x = start_x_dat, y = start_y_dat,
                     xend = end_x_dat, yend = end_y_dat),
                 arrow = arrow(length = unit(arrowlength, "cm")),
                 size = arrowsize, color = "black") +
    scale_color_gradient(low = "#f9ce34", high = "#6228d7") +
    labs(title = "Cell Pseudotime", color = "Pseudotime") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  return(list("ArrowOnPseudotime" = pseudotime_arrow_plot))
}

color_map <- c("#D57358","#8a508f","#ebe5c2", "#023047", "#F7CA71")

result <- plot_trajectory(pseudotime = Fig2E.data$pseudotime,
                          location =  cbind(-Fig2E.data$x, Fig2E.data$y),
                          clusterlabels =  Fig2E.data$subclone,
                          gridnum = 15, # adjust # of arrows
                          color_in = color_map,
                          pointsize = 0.3,
                          arrowdist = 50,
                          arrowlength = 0.15, # arrow tip length 
                          arrowsize = 0.45, # Arrow line thickness
                          textsize = 12)

print(result$ArrowOnPseudotime)

##----Figure 2f DEG----
library(Seurat)
library(tidyverse)
library(scRNAtoolVis)
library(RColorBrewer)

load('./Fig2f.RData')
p.vol <- jjVolcano(
  diffData = markers,
  tile.col = c("DCIS #1" = "#D57358", 
               "DCIS #2" = "#8a508f", 
               "Invasive" = "#023047", 
               "Other"   = "#ebe5c2"),
  size = 4.5,
  fontface = 'plain',
  aesCol = c("#1F78B4", "#E31A1C"),
  cluster.order = c("DCIS #1", "DCIS #2", "Invasive", "Other"),
  topGeneN = 5,
  polar = TRUE
)

p.vol <- p.vol +
  theme_void(base_size = 16) +  
  theme(
    plot.background = element_rect(fill = "white", color = NA),    
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),                                   
    axis.title = element_blank(),                           
    axis.text = element_blank(),                                 
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),                                
    legend.background = element_rect(fill = "white", color = "gray80")
  )

ggsave("./Figfg.png", plot = p.vol, width = 10, height = 10, dpi = 300)

pdf("./Fig2f.pdf", 10, 10)	 
print(p.vol)
dev.off()

##----Figure 2g Vocalno plot----
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggfun)
library(grid)
library(ggrepel)

### Invasive vs DCIS
load('./Fig2g1.RData')

p.i <- ggplot() +
  geom_point(
    data = gray_points_idc,
    aes(x = log2FoldChange, y = log10p_trunc, size = log10p_trunc),
    color = "gray70",
    alpha = 0.6
  ) +
  geom_point(
    data = color_points_idc,
    aes(x = log2FoldChange, y = log10p_trunc, size = log10p_trunc, color = log2FoldChange),
    alpha = 0.6
  ) +
  geom_point(
    data = top_genes_idc,
    aes(x = log2FoldChange, y = log10p_trunc, size = log10p_trunc, color = log2FoldChange),
    alpha = 1
  ) +
  geom_text_repel(
    data = top_genes_idc,
    aes(x = log2FoldChange, y = log10p_trunc, label = gene_id),
    size = 3.5,
    max.overlaps = 20,
    box.padding = 0.5,
    segment.ncp = 3,
    segment.angle = 20
  ) +
  scale_color_gradientn(
    colors = c("#39489f", "#39bbec", "#ffffbf", "#f38466", "#b81f25"),
    values = seq(0, 1, 0.2),
    limits = c(-5, 3),
    oob = scales::squish,
    name = expression(log[2]*"(Fold Change)")
  )+
  scale_size(
    range = c(2, 6),
    name = expression(-log[10]*"(p-value)")
  ) +
  scale_x_continuous(limits = c(-5, 3)) +
  scale_y_continuous(limits = c(0, 220), expand = expansion(mult = c(0, 0))) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.6) +
  geom_hline(yintercept = 100, lty = 4, col = "black", lwd = 0.6) +
  
  labs(
    title = "Invasive vs DCIS",
    x = expression(log[2]*"(Fold Change)"),
    y = expression(-log[10]*"(adjusted p-value)")
  ) +
  coord_cartesian(clip = "off") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.6, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave("./Fig2f_invasive_vs_dcis.png", plot = p.i,
       width = 8, height = 5, dpi = 300, units = "in")

pdf("./Fig2f_invasive_vs_dcis.pdf", 8, 5)	 
print (p.i)
dev.off()

### DCIS #2 vs DCIS #1
load('./Fig2g2.RData')

p.d <- ggplot() +
  
  geom_point(
    data = gray_points_idc,
    aes(x = log2FoldChange, y = log10p_trunc, size = log10p_trunc),
    color = "gray70",
    alpha = 0.6
  ) +
  
  geom_point(
    data = color_points_idc,
    aes(x = log2FoldChange, y = log10p_trunc, size = log10p_trunc, color = log2FoldChange),
    alpha = 0.6
  ) +

  geom_point(
    data = top_genes_idc,
    aes(x = log2FoldChange, y = log10p_trunc, size = log10p_trunc, color = log2FoldChange),
    alpha = 1
  ) +
  
  geom_text_repel(
    data = top_genes_idc,
    aes(x = log2FoldChange, y = log10p_trunc, label = gene_id),
    size = 3.5,
    max.overlaps = 20,
    box.padding = 0.5,
    segment.ncp = 3,
    segment.angle = 20
  ) +
  
  scale_color_gradientn(
    colors = c("#39489f", "#39bbec", "#ffffbf", "#f38466", "#b81f25"),
    values = seq(0, 1, 0.2),
    limits = c(-2, 3), 
    oob = scales::squish,
    name = expression(log[2]*"(Fold Change)")
  )+
  
  scale_size(
    range = c(2, 6),
    name = expression(-log[10]*"(p-value)")
  ) +
  
  scale_x_continuous(limits = c(-2, 3)) +
  scale_y_continuous(limits = c(0, 220), expand = expansion(mult = c(0, 0))) +

  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.6) +
  geom_hline(yintercept = 100, lty = 4, col = "black", lwd = 0.6) +
  
  labs(
    title = "DCIS #2 vs DCIS #1",
    x = expression(log[2]*"(Fold Change)"),
    y = expression(-log[10]*"(adjusted p-value)")
  ) +

  coord_cartesian(clip = "off") +
    theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.6, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave("volcano_DCIS2.png", plot = p.d,
       width = 8, height = 5, dpi = 300, units = "in")
ggsave("./Fig2f_dcis2_vs_dcis1.png", plot = p.d,
       width = 8, height = 5, dpi = 300, units = "in")

pdf("./Fig2f_dcis2_vs_dcis1.pdf", 8, 5)	 
print (p.i)
dev.off()

##----Figure 2h GSEA----
load('./Fig2h.RData')

p_list <- gsea.df %>%
  group_split(Clone) %>%
  lapply(function(data) {
    data <- data %>% arrange(desc(log10p))
    
    ggplot(data, aes(x = reorder(Term, log10p), y = log10p, fill = Clone)) +
      geom_col(width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("Invasive" = "#023047", "DCIS #1" = "#D57358", "DCIS #2" = "#8a508f")) +
      labs(x = NULL, y = expression(-log[10](p-value)), title = unique(data$Clone)) +
      theme(
        panel.background = element_blank(),       
        plot.background = element_blank(),        
        panel.grid.major = element_blank(),       
        panel.grid.minor = element_blank(),       
        panel.border = element_blank(),           
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = NA),  
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none"
      )
  })

final_plot <- p_list[[1]] + p_list[[2]] + p_list[[3]] + 
  plot_layout(ncol = 3)

ggsave("./Fig2h.png", final_plot, width = 15, height = 4, dpi = 300)

pdf("./Fig2h.pdf", 15, 4)	 
print(final_plot)
dev.off()

##----Figure 2i,j PGR/ESR1/ERBB2----
### Fig2i: z-score heatmap
load('./Fig2i.RData')

get_cluster_zscores <- function(seurat_obj, genes, group.by = "subclone") {
  
  clusters <- as.character(seurat_obj@meta.data[[group.by]])
  unique_clusters <- unique(clusters)
  
  mean_expr_mat <- matrix(
    nrow = length(genes),
    ncol = length(unique_clusters),
    dimnames = list(genes, unique_clusters)
  )
  
  expr_mat <- GetAssayData(seurat_obj, layer = "data")[genes, ]
  
  for (cluster in unique_clusters) {
    cells_in_cluster <- colnames(seurat_obj)[clusters == cluster]
    cluster_expr <- expr_mat[, cells_in_cluster, drop = FALSE]
    mean_expr_mat[, cluster] <- rowMeans(cluster_expr)
  }
  
  zscore_mat <- t(scale(t(mean_expr_mat)))
  
  return(list(
    mean_expression = mean_expr_mat,
    zscore = zscore_mat
  ))
}

results <- get_cluster_zscores(obj, genes_of_interest)
zscore_mat <- results$zscore
desired_order <- c("DCIS #1", "DCIS #2", "Invasive", "Other")
zscore_mat <- zscore_mat[, desired_order]
mean_mat <- results$mean_expression

pheatmap(
  mat = zscore_mat,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 30,
  cellheight = 30,
  fontsize = 7,
  main = "Z-score Expression",
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 7,
  filename = "./Fig2j.png",
  width = 10,
  height = 5
)

### Fig2j: expression heatmap
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

genes_of_interest <- c("PGR", "ESR1", "ERBB2")
obj$subclone <- factor(obj$subclone, levels = c("DCIS #1", "DCIS #2", "Invasive", "Other"))
Idents(obj) <- 'subclone'

subclone_colors <- c(
  "DCIS #1" = "#D57358",
  "DCIS #2" = "#8a508f",
  "Invasive" = "#023047",
  "Other"  = "#ebe5c2"
)

p <- DoHeatmap(
  object = obj,
  features = genes_of_interest,
  group.bar = TRUE,
  label = TRUE,
  size = 4,
  angle = 0,
  hjust = 0.5,
  draw.lines = TRUE,
  group.colors = subclone_colors
) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3")) +   
  scale_color_manual(                                                
    values = subclone_colors,
    name = "Cluster"                                                 
  ) +
  scale_y_discrete(limits = rev(genes_of_interest))                 

ggsave("./Fig2j.png", p, width = 10, height = 5.5, dpi = 300)

pdf("./Fig2j.pdf", 12, 4)	 
print(p)
dev.off()



