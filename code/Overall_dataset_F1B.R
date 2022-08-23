#R script to plot the overall distribution of the dataset
#across the sample sources and years

#clear the working environment

rm(list = ls())
graphics.off()

#load the require package(s)

require(ggplot2)
require(paletteer)

#select the colour palette for visualisation
colour <- paletteer_d("lisa::OskarSchlemmer")
#select the colour required
colour <- colour[c(1,3,5)]
#rearrange the colour
colour <- colour[c(3,2,1)]

#load in the dataset(s)
metadata <- read.csv("../data/BP.3341.metadata.source.edited2.csv", header = TRUE)

#ensure that the year is continuous variable, permitting the custom axis break 
metadata$year <- as.numeric(metadata$year)

breaks_custom <- seq(1935,2018, by = 5)

p <- ggplot(metadata, aes(x =as.factor(year), fill = Source)) +
  scale_fill_manual(name = "Source",
                    values = colour,
                    labels = c("Animal", "Environment", "Human", "NA"))+
  geom_bar(stat = "count") +
  scale_x_discrete(breaks = breaks_custom, labels = breaks_custom) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30,
                                   vjust = 1.5,
                                   hjust = 1),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 12.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "Year", 
       y = "Frequency") 

#save the plot
ggsave(plot = p, filename = "../results/BP_yr_distribution.pdf",
       width = 8.27, height = 5.83, units = "in" )
