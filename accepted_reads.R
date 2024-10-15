## create proportional bar chart showing number of accepted reads - post QC

library(ggplot2)

R_plot_QC_reads$Reads <- factor(R_plot_QC_reads$Reads, levels=c('Excluded_reads', 'Accepted_reads'))

myplot = ggplot(R_plot_QC_reads, aes(x = Samples, y= Counts, fill = Reads)) + 
  geom_col(colour = "black", position = "fill") +
  scale_fill_brewer(palette = "Set1") +
  geom_bar(position = 'stack', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, size = 10))

myplot + theme(panel.background = element_blank())
