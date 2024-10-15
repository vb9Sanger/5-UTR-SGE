library(ggplot2)

#create histogram displaying frequency of subsampled read lengths  
ggplot(histogram_lengths_data_all, aes(x = Read_Lengths)) + 
  geom_histogram(binwidth = 25, fill = "steelblue3", color = "black") + 
  
  labs(title = "HDR Library", x = "Read Lengths", y = "Frequency") + 
  theme_minimal()
