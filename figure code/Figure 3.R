# Download data here: https://drive.google.com/drive/folders/1DOcJTb7YFPwXm2Pa4tiZWPFzHxurqqLa?usp=drive_link

##Figure 3
#### Fig2b CNV heatmap: see tutorial website: https://czang409.github.io/highSpaClone/articles/xenium_breast2_code.html

##----Figure 3a Method Comparison----
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(scRNAtoolVis)

#6176 tumors
#358 bins

load('./F3a.RData')
make_cluster_plot <- function(data, colname, title, colors=c("#ebe5c2", "#D57358","#8a508f", "#023047", "#F7CA71")) {
  ggplot(data = data, aes(x = x, y = y, color = factor(.data[[colname]]))) +
    geom_point(size = 0.5) +
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

p1 <- make_cluster_plot(F3A.data, "highSpaClone", "highSpaClone")
p2 <- make_cluster_plot(F3A.data, "iris", "IRIS")
p3 <- make_cluster_plot(F3A.data, "infercnv", "InferCNV")
p4 <- make_cluster_plot(F3A.data, "truth", "Ground Truth")

# combine four figs
combined_plot <- p4 + p1 + p2 + p3 + plot_layout(ncol = 4)

# save figure
ggsave("./Fig3a.png", combined_plot, width = 20, height = 4, dpi = 300)
pdf("./Fig3a.pdf", 20, 4)
print (combined_plot)
dev.off()

##----Figure 2c CNV burden----
## total 358 bins
load('./F3c.RData')
custom_colors <- c("CNV_gain" = "#ffadad", "CNV_loss" = "#a0c4ff")
p <- ggplot(Fig3C.data, aes(x = highSpaClone, y = CNV_burdens, fill = CNV_status)) +
  geom_violin(position = position_dodge(width = 0.8),
              color = NA,
              alpha = 0.5,
              width = 0.8,
              trim = TRUE,
              scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6) +  # Remove outliers, reduce box width
  geom_point(data = Fig3C.data %>%
               group_by(highSpaClone, CNV_status) %>%
               sample_n(300),
             aes(x = highSpaClone, y = CNV_burdens, color = CNV_status),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             shape = 16,
             size = 1,
             alpha = 0.5,
             show.legend = FALSE)+

  scale_fill_manual(values = custom_colors) +  # Custom colors
  scale_x_discrete(labels = c("DCIS", "IDC")) +
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
  )+
  coord_cartesian(ylim = c(3, 25))

ggsave("./Fig3c.png", plot = p,
       width = 6.5, height = 4.5, dpi = 300, units = "in")

pdf("./Fig3c.pdf", 6.5, 4.5)
print(p)
dev.off()

##----Figure 3d trajectory----
load('./F3d.RData')

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
          axis.title = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal") +
    guides(color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm")))
  
  return(list("ArrowOnPseudotime" = pseudotime_arrow_plot))
}

result <- plot_trajectory(pseudotime = traj.df$pseudotime,
                          location =  cbind(traj.df$x, traj.df$y),
                          clusterlabels =  traj.df$subclone,
                          gridnum = 8, # adjust # of arrows
                          color_in = color_map,
                          pointsize = 0.3,
                          arrowdist = 50,
                          arrowlength = 0.15, # arrow tip length
                          arrowsize = 0.45, # Arrow line thickness
                          textsize = 12)

png("./Fig3d.png",width = 1800, height = 1200, res = 300)
print(result$ArrowOnPseudotime)
dev.off()

##----Figure 3e ERBB2----
load('./F3e.RData')
my_comparisons<- list( c("DCIS", "Other"),
                       c("IDC", "Other"),
                       c("IDC", "DCIS"))
p.vln <- VlnPlot(obj, features = 'ERBB2', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 4.5)) + NoLegend() +
  xlab(' ') + ggtitle('ERBB2') + ylab(' ') +
  scale_fill_manual(values = c('#d6e6ff', '#e5d4ef', '#ffadad', '#a0c4ff'))  +
  scale_x_discrete(labels = c('Other', 'IDC', 'DCIS'),
                   breaks = c('Other', 'IDC', 'DCIS'))
p.vln

ggsave("./Fig3e.png", plot = p.vln,
       width = 4, height = 4, dpi = 300, bg='white')

pdf("./Fig3e.pdf", 4, 4)
print(p)
dev.off()

##----Figure 3f Vocalno plot----
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggfun)
library(grid)
library(ggrepel)

###IDC vs DCIS
load('./F3f.RData')
p.idc <- ggplot() +
  
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
  
  scale_x_continuous(limits = c(-10, 3)) +
  scale_y_continuous(limits = c(0, 220), expand = expansion(mult = c(0, 0))) +
  
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.6) +
  geom_hline(yintercept = 100, lty = 4, col = "black", lwd = 0.6) +
  
  labs(
    title = "IDC vs DCIS",
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

ggsave("./Fig3f.png", plot = p.idc,
       width = 8, height = 5, dpi = 300, units = "in")

pdf("./Fig3f.pdf", 4, 4)
print(p.idc)
dev.off()

##----Figure 3g----
load('./F3g.RData')
set_colors <- function(pal, n) {
  
  if (all(pal %in% rownames(brewer.pal.info))) {
    num <- c()
    for (i in seq(length(pal))) {
      num[i] <- brewer.pal.info[pal[i],][[1]]
    }
    full_pal <- do.call(c, map2(.x = num, .y = pal, .f = brewer.pal))
  } else if (all(are_colors(pal))) {
    full_pal <- pal
  } else {
    stop('Incorrect palette setup. Please input valid RColorBrewer palette names or color names.')
  }
  
  if (n <= length(full_pal)) {
    return(full_pal[1:n])
  } else {
    warning("Number of colors required exceeds palette capacity. RdYlBu spectrum will be used instead.",
            immediate. = TRUE)
    return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
  }
}

plot_heatmap_gene_x_axis <- function(dataset,
                                     markers,
                                     sort_var = c('seurat_clusters'),
                                     n = 8,
                                     anno_var,
                                     anno_colors,
                                     gene_modules = NULL,  # gene module information
                                     module_colors = NULL, # module colors
                                     module_order = NULL,  # module display order
                                     gene_order_in_module = NULL,  # gene order within module
                                     sort_genes_by_module = TRUE, # whether to sort genes by module
                                     hm_limit = c(-2, 0, 2),
                                     hm_colors = c("#4575b4","white","#d73027"),
                                     row_font_size = 12) {
  
  require(ComplexHeatmap)
  require(circlize)
  require(dplyr)
  require(tibble)
  require(RColorBrewer)
  require(Seurat)
  
  # Get expression matrix
  mat <- GetAssayData(object = dataset, assay = DefaultAssay(dataset), slot = "scale.data")
  
  # Get genes
  if (is.data.frame(markers)) {
    genes <- get_top_genes(dataset, markers, n)
  } else if (is.character(markers)) {
    genes <- markers
  } else {
    stop('Incorrect input of markers')
  }
  
  genes <- intersect(genes, rownames(mat))
  
  # Sort genes by module
  if (!is.null(gene_modules) && sort_genes_by_module) {
    # Process gene module information
    if (is.data.frame(gene_modules)) {
      gene_module_vec <- setNames(gene_modules$module, gene_modules$gene)
    } else if (is.vector(gene_modules) && !is.null(names(gene_modules))) {
      gene_module_vec <- gene_modules
    } else {
      stop("gene_modules should be either a named vector or a data.frame with 'gene' and 'module' columns")
    }
    
    # Get module info for current genes
    current_gene_modules <- gene_module_vec[genes]
    current_gene_modules[is.na(current_gene_modules)] <- "Unknown"
    
    # Create sorting dataframe
    gene_order_df <- data.frame(
      gene = genes,
      module = current_gene_modules,
      stringsAsFactors = FALSE
    )
    
    # Set module order
    if (!is.null(module_order)) {
      # Use user-defined module order
      all_modules <- unique(current_gene_modules)
      missing_modules <- setdiff(all_modules, module_order)
      final_module_order <- c(module_order, missing_modules)  # Unspecified modules go last
      gene_order_df$module <- factor(gene_order_df$module, levels = final_module_order)
    } else {
      # Default: alphabetical order
      gene_order_df$module <- factor(gene_order_df$module)
    }
    
    # Sort genes within module by specified order
    if (!is.null(gene_order_in_module)) {
      # gene_order_in_module should be a named list, each element is a vector of gene order for one module
      # e.g. list("T_cell" = c("CD3D", "CD3E", "CD8A"), "B_cell" = c("MS4A1", "CD79A"))
      
      ordered_genes <- c()
      
      for (mod in levels(gene_order_df$module)) {
        # Get genes in current module
        genes_in_module <- gene_order_df$gene[gene_order_df$module == mod]
        
        if (mod %in% names(gene_order_in_module)) {
          # Use specified order
          specified_order <- gene_order_in_module[[mod]]
          # First genes in specified order, then unspecified genes (alphabetical)
          ordered_in_module <- intersect(specified_order, genes_in_module)
          unspecified_in_module <- sort(setdiff(genes_in_module, specified_order))
          ordered_genes <- c(ordered_genes, ordered_in_module, unspecified_in_module)
        } else {
          # No specified order, use alphabetical order
          ordered_genes <- c(ordered_genes, sort(genes_in_module))
        }
      }
      
      genes <- ordered_genes
    } else {
      # Sort by module, then by gene name (default logic)
      gene_order_df <- gene_order_df[order(gene_order_df$module, gene_order_df$gene), ]
      genes <- gene_order_df$gene
    }
  }
  
  mat <- mat[genes, , drop = FALSE]  # gene × cell
  
  # Get annotations, ensure matching order
  cells_use <- colnames(mat)
  anno <- dataset@meta.data %>%
    rownames_to_column(var = "barcode") %>%
    filter(barcode %in% cells_use) %>%
    arrange(match(barcode, cells_use))  # Ensure consistent order
  
  mat <- mat[, anno$barcode]  # Ensure column order matches
  
  # Transpose matrix: from gene × cell to cell × gene
  mat <- t(mat)  # Now cell × gene
  
  # Construct row annotations (cell annotations)
  row_annos <- list()
  for (i in seq_along(anno_var)) {
    value <- anno[[anno_var[i]]]
    err_msg <- paste('Incorrect specification for annotation colors for', anno_var[i])
    
    # Set factor levels order (if anno_colors is an ordered list)
    if (!is.numeric(value)) {
      if (is.list(anno_colors[[i]]) && !is.null(names(anno_colors[[i]]))) {
        # If colors are a named list, set factor levels by names
        value <- factor(value, levels = names(anno_colors[[i]]))
      } else {
        value <- factor(value)
      }
    }
    
    if (is.numeric(value)) {
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info)[brewer.pal.info$category != 'qual'])) {
        n_col <- brewer.pal.info[anno_colors[[i]], 'maxcolors'][[1]]
        pal <- brewer.pal(n = n_col, name = anno_colors[[i]])
        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)),
                              c(pal[2], pal[(n_col+1)%/%2], pal[n_col - 1]))
      } else if (length(anno_colors[[i]]) == 3 && all(are_colors(anno_colors[[i]]))) {
        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), anno_colors[[i]])
      } else {
        stop(err_msg)
      }
      
      ha <- rowAnnotation(a = value,
                          col = list(a = col_fun),
                          border = TRUE,
                          annotation_label = anno_var[i])
    } else {
      l <- levels(value)  
      
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info))) {
        col <- set_colors(anno_colors[[i]], length(l))
      } else if (length(anno_colors[[i]]) >= length(l) & all(are_colors(anno_colors[[i]]))) {
        col <- anno_colors[[i]]
      } else {
        stop(err_msg)
      }
      
      names(col) <- l
      col <- col[!is.na(names(col))]
      col <- list(a = col)
      
      ha <- rowAnnotation(a = value,
                          col = col,
                          border = TRUE,
                          annotation_label = "")
    }
    
    names(ha) <- anno_var[i]
    row_annos[[i]] <- ha
  }
  
  row_annos <- do.call(c, row_annos)
  row_annos@gap <- rep(unit(1, "mm"), length(row_annos))
  
  # Construct column annotations (gene module annotations)
  col_annos <- NULL
  if (!is.null(gene_modules)) {
    # gene_modules should be a named vector: names are gene names, values are module names
    # Or a data.frame with gene and module columns
    
    if (is.data.frame(gene_modules)) {
      # If data.frame format
      gene_module_vec <- setNames(gene_modules$module, gene_modules$gene)
    } else if (is.vector(gene_modules) && !is.null(names(gene_modules))) {
      # If named vector
      gene_module_vec <- gene_modules
    } else {
      stop("gene_modules should be either a named vector or a data.frame with 'gene' and 'module' columns")
    }
    
    # Get module info for current genes
    current_genes <- colnames(mat)
    gene_module_info <- gene_module_vec[current_genes]
    
    # Handle missing values
    gene_module_info[is.na(gene_module_info)] <- "Unknown"
    
    # Set module colors
    modules <- unique(gene_module_info)
    
    # Arrange modules by specified order
    if (!is.null(module_order)) {
      missing_modules <- setdiff(modules, module_order)
      modules <- c(intersect(module_order, modules), missing_modules)
    }
    
    if (is.null(module_colors)) {
      # Auto-generate colors
      if (length(modules) <= 8) {
        module_cols <- RColorBrewer::brewer.pal(min(8, length(modules)), "Set2")[1:length(modules)]
      } else {
        module_cols <- rainbow(length(modules))
      }
    } else {
      module_cols <- module_colors
    }
    
    names(module_cols) <- modules
    
    # create ordered factor to ensure legend order
    gene_module_info_ordered <- factor(gene_module_info, levels = modules)
    
    # Create column annotations
    col_annos <- columnAnnotation(
      Module = gene_module_info_ordered,  # ✅ Use ordered factor
      col = list(Module = module_cols),
      border = TRUE,
      #annotation_label = "",
      annotation_label = "Gene Module",
      annotation_name_gp = gpar(fontsize = 10)
    )
  }
  
  # Construct heatmap object
  ht <- Heatmap(mat,
                cluster_columns = FALSE,  # Do not cluster genes (columns)
                cluster_rows = TRUE,      # Cluster cells (rows) within groups
                row_split = factor(anno[[sort_var[1]]]),
                border = FALSE,
                rect_gp = gpar(col = NA),
                left_annotation = row_annos,
                top_annotation = col_annos,  # ✅ Add gene module annotations
                col = colorRamp2(hm_limit, hm_colors),
                heatmap_legend_param = list(direction = "horizontal",
                                            legend_width = unit(6, "cm"),
                                            title = "Expression"),
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = row_font_size))
  
  draw(ht,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "right")
}

png("./Fig3g.png",
    width = 8000, height = 6000, res = 600)

plot_heatmap_gene_x_axis(
  dataset = tumor.obj,
  markers = top20$gene,
  sort_var = "subclone",
  anno_var = "subclone",
  anno_colors =  list("Paired"),
  module_order = module_orders,
  module_colors = c( '#ffadad',  '#ffd6a5',  '#fdffb6',  '#caffbf',
                     '#9bf6ff',  '#a0c4ff',  '#bdb2ff',  '#ffc6ff',
                     '#008585', '#ba8274'),
  sort_genes_by_module = TRUE,
  gene_modules = module_df,
  gene_order_in_module = gene_order_in_module
)

dev.off()

##----Figure 3h GSEA----
load('./F3h.RData')
p <-  ggplot(gsea.df, aes(x = Term, y = log10p, fill = Clone)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("IDC" = "#8a508f")) +
  labs(x = NULL, y = expression(-log[10](p-value)), title = unique(gsea.df$Clone)) +
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

ggsave("./Fig3h.png", plot = p,
       width = 6, height = 4, dpi = 300, bg='white')

pdf("./Fig3h.pdf", 4, 4)
print(p)
dev.off()

##----Figure 3i----
load('./F3i.RData')
colors <- c("#D57358", "#023047", "#F7CA71","#8a508f","#ebe5c2")

p <- ggplot(data = data, aes(x = x, y = y, color = factor(pathology))) +
  geom_point(size = 0.5) +
  scale_color_manual(name = "Cluster", values = colors) +
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
    plot.margin       = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("./Fig3i1.png", plot = p,
       width = 8, height = 6, dpi = 300, units = "in")

pdf("./Fig3i1.pdf", 8, 6)
print(p)
dev.off()

####TSNE
load('./F3e.RData')
Idents(obj) <- 'pathology'

p <- DimPlot(obj, reduction = "tsne", pt.size = 1, cols=c("#ebe5c2","#8a508f", "#023047", "#F7CA71", "#D57358"))
ggsave("./Fig3k.png", plot = p,
       width = 8, height = 6, dpi = 300, units = "in")

ggsave("./Fig3i2.png", plot = p,
       width = 8, height = 6, dpi = 300, units = "in")

pdf("./Fig3i2.pdf", 8, 6)
print(p)
dev.off()

##----Figure 3j----
load('./F3j.RData')

p <- jjDotPlot(object = obj,
               gene = top_genes_unique,
               id = 'pathology',
               ytree = F)

ggsave("./Fig3j.png", plot = p,
       width = 12, height = 7.5, dpi = 300, units = "in")

pdf("./Fig3j.pdf", 12, 7.5)
print(p)
dev.off()

