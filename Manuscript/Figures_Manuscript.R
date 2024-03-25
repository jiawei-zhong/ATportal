setwd("~/05_Data/19_Portal/Figure_Prep")
counts <- read.delim("Count_Datasets.txt")
library(ggplot2)
counts$Module <- factor(counts$Module, levels = counts$Module)
ggplot(counts, aes(x=Module, y=Trancriptome, fill=Module))+ geom_col() + theme_bw() + scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = c("#4f736a", "#e2c744","#1a4659","#ba5c28", "#f59b7c", "#b2dfda"))
ggplot(counts, aes(x=Module, y=Total, fill=Module))+ geom_col() + theme_bw() + scale_y_continuous(expand=c(0,0)) + scale_x_discrete(limits=rev)+
  scale_fill_manual(values = c("#4f736a", "#e2c744","#1a4659","#ba5c28", "#f59b7c", "#b2dfda")) + coord_flip() 
