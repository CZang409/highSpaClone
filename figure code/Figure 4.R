# Download data here: https://drive.google.com/drive/folders/1EG_I6DjPyuD2bIXXGPR7OJwFsfpB5I9w?usp=drive_link

##Figure 4
####Fig4c and 4e: see tutorial website: https://czang409.github.io/highSpaClone/articles/Intro_to_highSpaClone.html
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)

##----Figure 4a Cell type annotation----
colors <- c('#00202e',  '#003f5c',  '#2c4875',  '#8a508f',  '#bc5090',  '#ff6361',  '#ff8531',  '#ffa600', '#ffd380', '#78938a')

load('./F4.RData')
p <- ggplot(data = data, aes(x = x, y = y, color = factor(celltype))) +
  geom_point(size = 0.01) +
  scale_x_reverse()+
  scale_color_manual(name = "Cell type", values = colors)+ 
  coord_flip()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              face = 'bold',
                              size = 16),
    plot.title.position = "plot",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.border = element_rect(color = "black", fill = NA),
    plot.background = element_rect(fill = NA),
  ) +
  ggtitle("")+
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("./F4a.png", p, width = 5, height = 4.5, dpi = 300)
pdf("./Fig4a.pdf", 5, 4.5)	 
print(p)
dev.off()

##----Figure 4b Tumor cell annotation----
load('./F4.RData')
colors <- c("#ebe5c2", "#D57358")

p <- ggplot(data = data, aes(x = x, y = y, color = factor(tumor))) +
  geom_point(size = 0.05) +
  scale_x_reverse()+
  scale_color_manual(name = "Cell Type", values = colors)+ 
  coord_flip()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              face = 'bold',
                              size = 16),
    plot.title.position = "plot",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.border = element_rect(color = "black", fill = NA),
    plot.background = element_rect(fill = NA),
  ) +
  ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("./F4b.png", p, width = 5, height = 4.5, dpi = 300)
pdf("./Fig4b.pdf", 5, 4.5)	 
print(p)
dev.off()

##----Figure 4d Tumor subclone visualization----
load('./F4.RData')
colors <- c("#E64B35", "#4DBBD5", "#8491B4","#ebe5c2")

p <- ggplot(data = data, aes(x = x, y = y, color = factor(subclone))) +
  geom_point(size = 0.05) +
  scale_x_reverse()+
  scale_color_manual(name = "Cluster", values = colors)+ 
  coord_flip()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              face = 'bold',
                              size = 16),
    plot.title.position = "plot",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.border = element_rect(color = "black", fill = NA),
    plot.background = element_rect(fill = NA),
  ) +
  ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("./F4d.png", p, width = 5, height = 4.5, dpi = 300)
pdf("./Fig4d.pdf", 5, 4.5)	 
print(p)
dev.off()


##----Figure 4f CNV burden boxplot----
load('./F4f.RData')
custom_colors <- c("CNV_gain" = "#ffadad", "CNV_loss" = "#a0c4ff")

cnv_filt <- cnv_burden_long_data %>%
  group_by(cell.label, CNV_status) %>%
  mutate(
    Q1 = quantile(CNV_burdens, 0.25, na.rm = TRUE),
    Q3 = quantile(CNV_burdens, 0.75, na.rm = TRUE),
    IQRv = Q3 - Q1,
    lo = Q1 - 1.5 * IQRv,
    hi = Q3 + 1.5 * IQRv
  ) %>%
  ungroup() %>%
  filter(CNV_burdens >= lo, CNV_burdens <= hi)

p <- ggplot(cnv_burden_long_data, aes(x = cell.label, y = CNV_burdens, fill = CNV_status)) +
  geom_violin(position = position_dodge(width = 0.8),
              color = NA,
              alpha = 0.5,
              width = 0.8,
              trim = TRUE,
              scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6) +  # Remove outliers, reduce box width
  geom_point(data = cnv_burden_long_data %>%
               group_by(cell.label, CNV_status) %>%
               sample_n(300),
             aes(x = cell.label, y = CNV_burdens, color = CNV_status),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             shape = 16,
             size = 1,
             alpha = 0.5,
             show.legend = FALSE)+
  
  scale_fill_manual(values = custom_colors) +  # Custom colors
  scale_x_discrete(labels = c("Clone 1", "Clone 2", "Clone 3")) +
  labs(title = "CNV Burdens Across Clones",
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
  coord_cartesian(ylim = c(10, 50))
p

ggsave("./F4f.png", p, width = 6.5, height = 4.5, dpi = 300)
pdf("./Fig4f.pdf", 6.5, 4.5)	 
print(p)
dev.off()

##----Figure 4g ESR1, ERBB2, PGR markers----
load('./F4g.RData')

my_comparisons <- list(
  c("Clone 1", "Clone 2"),
  c("Clone 1", "Clone 3"),
  c("Clone 2", "Clone 3")
)

p.vln <- VlnPlot(obj_filtered2, features = 'ERBB2', pt.size = 0,
                 cols = c('#d6e6ff', '#e5d4ef', '#ffadad')) +  
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2, fill = NA) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label = "p.signif", method = "wilcox.test") +
  coord_cartesian(ylim = c(0, 7)) + 
  NoLegend() +
  labs(title = "ERBB2", x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

ggsave("./ERBB2_violin.png", p.vln, width = 4.5, height = 3, dpi = 300, bg = 'white')
pdf("./ERBB2_violin.pdf", 4.5, 3)	 
print(p.vln)
dev.off()

p.vln <- VlnPlot(obj_filtered2, features = 'ESR1', pt.size = 0,
                 cols = c('#d6e6ff', '#e5d4ef', '#ffadad')) +  
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2, fill = NA) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label = "p.signif", method = "wilcox.test") +
  coord_cartesian(ylim = c(0, 4)) + 
  NoLegend() +
  labs(title = "ESR1", x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

ggsave("./ESR1_violin.png", p.vln, width = 4.5, height = 3, dpi = 300, bg = 'white')
pdf("./ESR1_violin.pdf", 4.5, 3)	 
print(p.vln)
dev.off()


p.vln <- VlnPlot(obj_filtered2, features = 'PGR', pt.size = 0,
                 cols = c('#d6e6ff', '#e5d4ef', '#ffadad')) +  
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2, fill = NA) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label = "p.signif", method = "wilcox.test") +
  coord_cartesian(ylim = c(0, 4)) + 
  NoLegend() +
  labs(title = "PGR", x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

ggsave("./PGR_violin.png", p.vln, width = 4.5, height = 3, dpi = 300, bg = 'white')
pdf("./PGR_violin.pdf", 4.5, 3)	 
print(p.vln)
dev.off()


##----Figure 4h TLS score----
# calculate tls score
tls.markers <- list(c('IGHA1', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 
                      'IGHGP', 'IGHM', 'IGKC', 'IGLC1', 'IGLC2', 
                      'IGLC3', 'JCHAIN', 'CD79A', 'FCRL5', 'MZB1', 'SSR4', 
                      'XBP1', 'TRBC2', 'IL7R', 'CXCL12', 'LUM', 'C1QA', 'C7', 
                      'CD52', 'APOE', 'PLTP', 'PTGDS', 'PIM2', 'DERL3'))
obj <- AddModuleScore(
  object = obj,
  features = tls.markers,
  name = 'TLS.score'
)

# visualization
load('./F4.RData')

p <- ggplot() +
  geom_point(
    data = subset(data, tumor == "Tumor"),
    aes(x = x, y = y),
    color = "gray80",
    alpha = 0.9,
    size = 0.05
  ) +
  geom_point(
    data = subset(data, tumor != "Tumor"),
    aes(x = x, y = y, color = tls.score),
    size = 0.05
  ) +
  scale_x_reverse() +
  scale_color_gradient(
    name = "TLS.score",
    low = "#6228d7",
    high = "#f9ce34",
    limits = c(0, 2),
    breaks = seq(0, 2, by = 0.5),
  ) +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              face = 'bold',
                              size = 16),
    plot.title.position = "plot",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.border = element_rect(color = "black", fill = NA),
    plot.background = element_rect(fill = NA),
    legend.position = "bottom",                  
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold")
  ) +
  ggtitle("") 

ggsave("./F4h.png", p, width = 3.5, height = 4.5, dpi = 300)
pdf("./Fig4h.pdf", 3.5, 4.5)	 
print(p)
dev.off()


##----Figure 4i DEG----
load('./F4i.RData')

p <- ggplot(data, aes(x = features.plot, y = id)) +
  geom_tile(fill = NA, color = "black", linewidth = 0.3) +
  geom_point(aes(size = pct.exp, color = expr01)) +
  scale_color_gradientn(colours = c("#0571b0", "#92c5de", "#f4a582", "#ca0020"),
                        limits = c(0, 1), name = "Average Expression") +
  scale_size(range = c(1, 6), limits = c(0, 100), name = "Percent Expressed",
             breaks = c(25, 50, 75, 100)) +
  
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    legend.position = "top",          
    legend.box = "horizontal",        
    legend.direction = "horizontal",  
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(1, "cm"), 
    plot.margin = margin(5, 5, 5, 5, "mm")
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",         
      title.hjust = 0.5,             
      barwidth = unit(3, "cm"),       
      barheight = unit(0.35, "cm"),    
      ticks.colour = "black",
      ticks.linewidth = 0.3
    ),
    size = guide_legend(
      title.position = "top",         
      title.hjust = 0.5,          
      nrow = 1,                      
      override.aes = list(color = "grey30")
    )
  )

ggsave("./F4i.png", p, width = 6.5, height = 2.5, dpi = 300)
pdf("./Fig4d.pdf", 6.5, 2.5)	 
print(p)
dev.off()

##----Figure 4j GSEA----
load('./F4j.RData')

p <- ggplot(data, aes(x = Clone, y = Pathway, fill = NES)) +
  geom_tile(color = "#8ba0a4") +
  scale_fill_gradient2(low = "#b8e0d4", mid = "white", high = "#f27a7d", midpoint = 0, na.value = '#d1dbe4') +
  labs(title = "", x = "", y = "", fill = "NES") +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1))

ggsave("./F4j.png", p, width = 8, height = 4, dpi = 300)
pdf("./Fig4j.pdf", 8, 4)	 
print(p)
dev.off()

