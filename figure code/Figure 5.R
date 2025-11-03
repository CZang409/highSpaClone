# Download data here: https://drive.google.com/drive/folders/1-djdJ_SBWkTA9p3bWWM7hTUE1Q3AzGbe?usp=drive_link 

##Figure 5
####Fig5c CNV heatmap: see tutorial website: https://czang409.github.io/highSpaClone/articles/visiumhd_colon_code.html
library(ggplot2)
library(dplyr)
library(patchwork)
library(VennDiagram)
library(tidyr)

##----Figure 5a Cell type composition----
load('./F5a.RData')
make_cluster_plot <- function(data, colname, title='', colors=c('#00202e',  '#003f5c',  '#2c4875',  '#8a508f',  '#bc5090',  '#ff6361',  '#ff8531',  '#ffa600', '#ffd380', '#78938a')) {
  ggplot(data = data, aes(x = x, y = y, color = factor(.data[[colname]]))) +
    geom_point(size = 0.01) +
    scale_x_reverse() +
    coord_flip()+
    scale_color_manual(name = "Cell type", values = colors) +
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

p1 <- make_cluster_plot(label1, "celltype", 'P1CRC')
p2 <- make_cluster_plot(label2, "celltype", 'P2CRC')
p3 <- make_cluster_plot(label3, "celltype", 'P5CRC')
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave("./Fig5a.png", combined_plot, width = 12, height = 4, dpi = 300)

pdf("./Fig5a.pdf", 20, 4)	 
print (combined_plot)
dev.off()

##----Figure 5b Tumor visualization----
load('./F5a.RData')
colors <- c("#8a508f","#ebe5c2", "#D57358")

p1 <- make_cluster_plot(label1, "tumor", 'P1CRC',colors=colors)
p2 <- make_cluster_plot(label2, "tumor", 'P2CRC',colors=colors)
p3 <- make_cluster_plot(label3, "tumor", 'P5CRC',colors=colors)
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave("./Fig5b_new.png", combined_plot, width = 12, height = 4, dpi = 300)

pdf("./Fig5b.pdf", 20, 4)	 
print (combined_plot)
dev.off()

##----Figure 5d CNV score----
load('./F5d.RData')

plot_cnv_score <- function(data, 
                           colname,
                           title='') {
  ggplot(data = data, aes(x = x, y = y, color = .data[[colname]])) +
    geom_point(size = 0.01) +
    scale_x_reverse() +
    scale_color_gradient(
      low = "#e5d4ef", 
      high = "#8a508f", 
      limits = c(0.35, 0.5),
      na.value = "gray80"
    ) +
    labs(title = title, color = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "bottom"
    ) +
    coord_flip() +
    guides(color = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(3, "cm"),
      barheight = unit(0.4, "cm"),
      title.position = "top",
      title.hjust = 0.5
    ))
}

p1 <- plot_cnv_score(p1.cnv, 'cnv.score.capped', title = 'P1CRC')
p2 <- plot_cnv_score(p2.cnv, 'cnv.score.capped', title = 'P2CRC')
p3 <- plot_cnv_score(p5.cnv, 'cnv.score.capped', title = 'P5CRC')
combined_plot <- p1 + p2 + p3

ggsave("./Fig5d.png", combined_plot, width = 12, height = 4.5, dpi = 300)

pdf("./Fig5d.pdf", 12, 4.5)	 
print (combined_plot)
dev.off()

##----Figure 5d venn plot----
load('./F5e.RData')

set1 <- names(cnv_status.p1[cnv_status.p1 != 0])
set2 <- names(cnv_status.p2[cnv_status.p2 != 0])
set3 <- names(cnv_status.p5[cnv_status.p5 != 0])

venn.plot <- venn.diagram(
  x = list(P1CRC = set1, P2CRC = set2, P5CRC = set3),
  filename = NULL,
  
  main = "Shared and Unique CNV Bins in Tumor Section",
  main.cex = 1.6,
  
  lwd = 1.2,
  col = "#4D4D4D",
  fill = c("#66c2a5", "#fc8d62", "#8da0cb"),
  alpha = 0.55,          
  margin = 0.06,           
  
  cex = 1.4,
  label.col = "#333333",
  fontfamily = "Helvetica",

  category.names = c("P1CRC", "P2CRC", "P5CRC"),
  cat.cex = 1.6,
  cat.fontface = "bold",
  cat.col = "black",
  cat.default.pos = "outer",
  cat.pos  = c(-15, 15, 0), 
  cat.dist = c(0.10, 0.10, 0.10),
  cat.just = list(c(1,0.5), c(0,0.5), c(0.5,0.5)),
  
  output = TRUE,
  imagetype = "png",
  height = 2000, width = 2000, resolution = 500
)

grid.newpage()
grid.draw(venn.plot)

##----Figure 5f CNV burden----
load('./F5f.RData')
## transite wide data to long data
cnv_long_data <- cnv.burden %>%
  pivot_longer(cols = c("CNV_gain", "CNV_loss"), 
               names_to = "CNV_status", 
               values_to = "CNV_burdens")

custom_colors <- c("CNV_gain" = "#ffadad", "CNV_loss" = "#a0c4ff")

cnv_filt <- cnv_long_data %>%
  group_by(sample, CNV_status) %>%
  mutate(
    Q1 = quantile(CNV_burdens, 0.25, na.rm = TRUE),
    Q3 = quantile(CNV_burdens, 0.75, na.rm = TRUE),
    IQRv = Q3 - Q1,
    lo = Q1 - 1.5 * IQRv,
    hi = Q3 + 1.5 * IQRv
  ) %>%
  ungroup() %>%
  filter(CNV_burdens >= lo, CNV_burdens <= hi)

p <- ggplot(cnv_filt, aes(x = sample, y = CNV_burdens, fill = CNV_status)) +
  geom_violin(position = position_dodge(width = 0.8),
              color = NA,
              alpha = 0.5,
              width = 0.8,
              trim = TRUE,
              scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6) +  # Remove outliers, reduce box width
  geom_point(data = cnv_filt %>%
               group_by(sample, CNV_status) %>%
               sample_n(300),
             aes(x = sample, y = CNV_burdens, color = CNV_status),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             shape = 16,
             size = 1,
             alpha = 0.5,
             show.legend = FALSE)+
  
  scale_fill_manual(values = custom_colors) +  # Custom colors
  scale_x_discrete(labels = c("P1CRC", "P2CRC", "P5CRC")) +
  labs(title = "CNV Burdens Across Patients",
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
  coord_cartesian(ylim = c(0, 90))

ggsave("./Fig5f.png", plot = p,
       width = 6.5, height = 4.5, dpi = 300, units = "in")

pdf("./Fig5f.pdf", 6.5, 4.5)	 
print (p)
dev.off()


