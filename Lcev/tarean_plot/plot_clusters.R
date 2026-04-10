#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(viridis)
library(hash)

args <- commandArgs(trailingOnly = TRUE)
input_clusters <- args[1]
input_stats <- args[2]
output_clusters <- args[3]
output_stats <- args[4]

df <- read.table(input_clusters, header = TRUE, sep = "\t")
df2 <- read.table(input_stats, header = TRUE, sep = "\t")

df_sum <- df[order(-df$Number.of.reads), ]

ann_colors <- c(
  "LTR" = "#440154",
  "LINE" = "#ff10e3",
  "SINE" = "#21918C",
  "Helitron" = "#5DC863",
  "Mutator" = "#277DA1",
  "hAT" = "#F8961E",
  "TcMar" = "#F3722C",
  "PiggyBac" = "#F94144",
  "PIF-Harbinger" = "#9B5DE5",
  "CMC" = "#577590",
  "P" = "#43AA8B",
  "MITE" = "#4D908E",
  "TIR" = "#90BE6D",
  "Organelle" = "#FDE725",
  "rDNA" = "#B56576",
  "NA" = "#898989"
)

df_sum$Annotation[df_sum$Annotation == "" |
                     is.na(df_sum$Annotation) |
                     df_sum$Annotation == "unknown"] <- "NA"

df_sum$Annotation[is.na(df_sum$Annotation)| df_sum$Annotation=="<NA>"] <- "NA"

df_sum$Cluster <- factor(df_sum$Cluster, levels = df_sum$Cluster)

df_sum$Annotation <- factor(df_sum$Annotation, levels = names(ann_colors))
df_sum$Annotation[is.na(df_sum$Annotation)] <- "NA"

p <- ggplot(df_sum, aes(
  x = Cluster,
  y = Number.of.reads,
  fill = Annotation
)) +
  geom_col() +
  labs(
    x = "Clusters",
    y = "Cluster size (n° of reads)"
  ) +
  theme(
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold", color = "black")
  ) +
  scale_fill_manual(values = ann_colors, drop = TRUE)

q <- ggplot(df2, aes(x = "", y = Proportion, fill = Type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
      geom_text(
        aes(label = ifelse(Proportion > 1.0,
                   sprintf("%.1f", Proportion),
                   "")),
        position = position_stack(vjust = 0.5),
        color = "black", size = 5,
        fontface = "bold") +
    labs(labs(title = "", x = "", y = "")) +
    theme(
        panel.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold", size = 15),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
    )+
  scale_fill_manual(values = ann_colors, drop = TRUE)

ggsave(output_stats, q, width = 8, height = 6, dpi = 300)
