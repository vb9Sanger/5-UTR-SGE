## create genomic coverage plots for mapped reads 

library(ggplot2)

########## Create a scatterplot displaying Log2(count+1) versus position for all variants 

#prepare data 
FINAL_count_positions_LOG2$VALUE <- as.numeric(FINAL_count_positions_LOG2$VALUE)

FINAL_count_positions_LOG2$POSITION<- as.numeric(FINAL_count_positions_LOG2$POSITION)

###create the plot 
plot1 <- ggplot(FINAL_count_positions_LOG2, aes(x=POSITION, y=VALUE)) + geom_point() + xlim(0,200) + ylim(0,20) + geom_point(color='mediumblue') + xlab("variant position") + ylab("log2(count+1)") + ggtitle("HDR Library")

plot1 + geom_hline(yintercept=2.59, linetype="dashed", 
                   color = "red3",)

