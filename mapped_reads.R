## create proportional bar chart showing % mapped reads - including WT and PAM only 

R_plot_QC_reads_WT$Reads <- factor(R_plot_QC_reads_WT$Reads, levels=c('Unmapped_reads', 'WT_reads', 'PAM_only_reads', 'Library_reads'))

myplot = ggplot(R_plot_QC_reads_WT, aes(x = Samples, y= Percent, fill = Reads)) + 
  geom_col(colour = "black", position = "fill") +
  scale_fill_manual(values = c("red3", "seagreen3", "sienna2", "steelblue3")) +
  geom_bar(position = 'stack', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, size = 10))

myplot + theme(panel.background = element_blank())