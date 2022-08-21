#Clean the working environment
rm(list = ls())
graphics.off()

#Load require package(s)
require(ggplot2)
require(viridisLite)

#Load the dataset(s)
metadata <- read.csv("../data/BP.3341.metadata.source.edited2.csv", header = TRUE)

metadata$year <- as.numeric(metadata$year)

#create a custom x axis breaks
year_break <- seq(min(metadata$year, na.rm = TRUE), max(metadata$year, na.rm = TRUE), by = 5)

p <- ggplot(metadata, aes(x =year, fill = Source)) +
  scale_fill_viridis_d(option = "D", na.value = "grey50", direction = -1,
                       labels = c("Animal", "Environment", "Human", "NA"))+
  geom_bar(stat = "count") +
  scale_x_continuous(breaks = year_break, labels = year_break) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30,
                                   vjust = 1.5,
                                   hjust = 1),
        legend.text = element_text(size = 12.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "Year", 
       y = "Frequency") 

#save the plot
ggsave(plot = p, filename = "../results/BP_yr_distribution.pdf",
       width = 8.27, height = 5.83, units = "in" )