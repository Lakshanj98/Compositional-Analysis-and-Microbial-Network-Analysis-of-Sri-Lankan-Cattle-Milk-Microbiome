library(ggplot2)

setwd("D:/Research/R codes/Functional Annotations/otu abundance approach/network functions")

###################################################
# generate graphs for entire networks_ function and different OTU count

# read the data
full_net <- read.table("full net network functions and otu abundance.tsv", sep = "\t", header = TRUE)
cl1_net <- read.table("high c network functions and otu abundance.tsv", sep = "\t", header = TRUE)
cl2_net <- read.table("medium c network functions and otu abundance.tsv", sep = "\t", header = TRUE)
cl3_net <- read.table("low c network functions and otu abundance.tsv", sep = "\t", header = TRUE)
large_net <- read.table("large fs network functions and otu abundance.tsv", sep = "\t", header = TRUE)
med_net <- read.table("med small network functions and otu abundance.tsv", sep = "\t", header = TRUE)

########## -------------------------------------------------------------------------
########## ------------- Create the bar plot for full network
full_plot <-ggplot(full_net, aes(x = OTU_abundance, y = Function)) +
  geom_bar(stat = "identity" , fill = "#00b300") +
  labs(title = "Functional frequencies in full network",
       x = "Functional Frequency",
       y = "Function") +
  theme_classic() 
#scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 100)) 

plot(full_plot)
# Adjust the plot appearance (optional)
theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))  # Adjust margins for landscape
# Save the plot (optional)
ggsave("Full net functional frequencies.pdf", width = 10, height = 7, units = "in")  # Adjust dimensions and format as needed

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Full net functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution

########## -------------------------------------------------------------------------
########## ------------- Create the bar plot for cl1_75_95 network
full_plot <-ggplot(cl1_net, aes(x = OTU_abundance, y = Function)) +
  geom_bar(stat = "identity" , fill = "#00b300") +
  labs(title = "Functional frequencies in high cleanliness network",
       x = "Functional Frequency",
       y = "Function") +
  theme_classic() 
#scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 100)) 

plot(full_plot)
# Adjust the plot appearance (optional)
theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))  # Adjust margins for landscape
# Save the plot (optional)
ggsave("High c net functional frequencies.pdf", width = 10, height = 6, units = "in")  # Adjust dimensions and format as needed

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("High c net functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution

########## -------------------------------------------------------------------------
########## ------------- Create the bar plot for cl2_55_75 network
full_plot <-ggplot(cl2_net, aes(x = OTU_abundance, y = Function)) +
  geom_bar(stat = "identity" , fill = "#00b300") +
  labs(title = "Functional frequencies in medium cleanliness network",
       x = "Functional Frequency",
       y = "Function") +
  theme_classic() 
#scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 100)) 

plot(full_plot)
# Adjust the plot appearance (optional)
theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))  # Adjust margins for landscape
# Save the plot (optional)
ggsave("Medium c net functional frequencies.pdf", width = 10, height = 7, units = "in")  # Adjust dimensions and format as needed

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Medium c net functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution

########## -------------------------------------------------------------------------
########## ------------- Create the bar plot for cl3_35_55 network
full_plot <-ggplot(cl3_net, aes(x = OTU_abundance, y = Function)) +
  geom_bar(stat = "identity" , fill = "#00b300") +
  labs(title = "Functional frequencies in low cleanliness network",
       x = "Functional Frequency",
       y = "Function") +
  theme_classic() 
#scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 100)) 

plot(full_plot)
# Adjust the plot appearance (optional)
theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))  # Adjust margins for landscape
# Save the plot (optional)
ggsave("Low c net functional frequencies.pdf", width = 10, height = 7, units = "in")  # Adjust dimensions and format as needed

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Low c functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution

########## -------------------------------------------------------------------------
########## ------------- Create the bar plot for large farm size network
full_plot <-ggplot(large_net, aes(x = OTU_abundance, y = Function)) +
  geom_bar(stat = "identity" , fill = "#00b300") +
  labs(title = "Functional frequencies in large farm size network",
       x = "Functional Frequency",
       y = "Function") +
  theme_classic() 
#scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 100)) 

plot(full_plot)
# Adjust the plot appearance (optional)
theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))  # Adjust margins for landscape
# Save the plot (optional)
ggsave("Large net functional frequencies.pdf", width = 10, height = 7, units = "in")  # Adjust dimensions and format as needed

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Large net functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution

########## -------------------------------------------------------------------------
########## ------------- Create the bar plot for med small farm size network
full_plot <-ggplot(med_net, aes(x = OTU_abundance, y = Function)) +
  geom_bar(stat = "identity" , fill = "#00b300") +
  labs(title = "Functional frequencies in med_small farm size network",
       x = "Functional Frequency",
       y = "Function") +
  theme_classic() 
#scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 100)) 

plot(full_plot)
# Adjust the plot appearance (optional)
theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))  # Adjust margins for landscape
# Save the plot (optional)
ggsave("Med_small net functional frequencies.pdf", width = 10, height = 7, units = "in")  # Adjust dimensions and format as needed

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Med_small net functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution

########## -------------------------------------------------------------------------
