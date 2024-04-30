# Load required library
library(ggplot2)

setwd("D:/Research/R codes/Functional Annotations/otu abundance approach/network cluster functions")

colors <- c(
  "#402010", "#44FFFF", "#FF0000", "#00FF00", "#0000FF",  
  "#FFFF00", "#FFA500", "#8B0080", "#FFC0CB", "#808080",  
  "#008080", "#000080", "#800000", "#808000", "#C000C0",  
  "#FF00FF", "#00FF00", "#00FFFF", "#4169E1", "#FF7F50",  
  "#ADD8E6", "#FFC0CB", "#90EE90", "#FFFFE0", "#E6E6FA",  
  "#87CEEB", "#98FB98", "#FFE5B4", "#87CEEB", "#00E9F5",  
  "#FFFFDD", "#000040", "#af0000", "#C0C0C0", "#BC8F8F",
  "#C39797", "#DDA0DD", "#BDB7B5", "#8B4513", "#FFA07A"
)

colors <- c(
  "#402010", "#FFFFFF", "#FF0000", "#00FF00", "#0000FF",  
  "#FFFF00", "#FFA500", "#8B0080", "#FFC0CB", "#808080",  
  "#008080", "#000080", "#800000", "#808000", "#C000C0",  
  "#FF00FF", "#00FF00", "#00FFFF", "#4169E1", "#FF7F50",  
  "#ADD8E6", "#FFC0CB", "#90EE90", "#FFFFE0", "#E6E6FA",  
  "#87CEEB", "#98FB98", "#FFE5B4", "#87CEEB", "#00E9F5",  
  "#FFFFDD", "#000040", "#FFD700", "#C0C0C0", "#BC8F8F",
  "#C39797", "#DDA0DD", "#BDB7B5", "#8B4513", "#FFA07A"
)
# used color codes
colors <- c(
  "#FF0000", "#FF8000", "#FFFF00", "#00FF00", "#0000FF",  
  "#00FFFF", "#800080", "#808080", "#C0C0C0", "#FF00FF",  
  "#AAFFFF", "#000080", "#800000", "#FFFFE0", "#F0FFF0",  
  "#008080", "#C39797", "#FFC0CB", "#DDA0DD", "#BDB7B5",  
  "#E6E6FA", "#FFF0F5", "#CD5B4B", "#F0E68C", "#EEE8AA",  
  "#98FB98", "#ADFF2F", "#FFD700", "#FFA500", "#FFFF99",  
  "#BC8F8F", "#4169E1", "#8B4513", "#FA8072", "#FFA07A",  
  "#A9A9A9", "#D3D3D3", "#F7F7F7", "#CCCCCC", "#F5F5F5"   
)

############################################################
#######       read all the data files
full_data <- read.table("full network cluster functions and otu abundance.tsv", sep = "\t", header = TRUE)

cl1_data<- read.table("cl1 network cluster functions and otu abundance.tsv", sep = "\t", header = TRUE)
cl2_data <- read.table("cl2 network cluster functions and otu abundance.tsv", sep = "\t", header = TRUE)
cl3_data <- read.table("cl3 network cluster functions and otu abundance.tsv", sep = "\t", header = TRUE)
large_data <- read.table("Large network cluster functions and otu abundance.tsv", sep = "\t", header = TRUE)
med_small_data <- read.table("Med network cluster functions and otu abundance.tsv", sep = "\t", header = TRUE)

########################################################### full network
# Create the bar plot with dodge placement and fill by function

my_plot = ggplot(full_data, aes(x = Cluster, y = OTU_abundance, fill = Function)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Functional frequencies in full network", x = "Cluster", y = "Functional Frequency") +
  theme_classic() +
  #scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 1))
  # Set fill color manually based on the unique function values
  scale_fill_manual(values = colors[match(full_data$Function, unique(full_data$Function))])+
  theme(legend.position = "right",           # Move legend to the right side
        legend.title = element_text(face = "bold"),  # Bold title
        legend.text = element_text(size = 7),legend.key.size = unit(0.3, "cm")) 

plot(my_plot)
theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Full network_clusters_Functional frequencies.pdf", width = 11, height = 8.5, units = "in")  # Set maximum paper size (adjust units as needed)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  # Save as JPEG with maximum paper size (adjust units as needed)
  ggsave("Full network_clusters_Functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution


############################################ cl1_75_95

# Create the bar plot with dodge placement and fill by function
my_plot = ggplot(cl1_data, aes(x = Cluster, y = OTU_abundance, fill = Function)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Functional frequencies in high cleanliness network", x = "Cluster", y = "Functional Frequency") +
  theme_classic() +
  #scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 1))
  scale_fill_manual(values = colors[match(cl1_data$Function, unique(cl1_data$Function))])+
  theme(legend.position = "right",           # Move legend to the right side
        legend.title = element_text(face = "bold"),  # Bold title
        legend.text = element_text(size = 7),legend.key.size = unit(0.3, "cm")) 
plot(my_plot)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("High clean net_clusters_Functional frequencies.pdf", width = 11, height = 8.5, units = "in")  # Set maximum paper size (adjust units as needed)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  # Save as JPEG with maximum paper size (adjust units as needed)
  ggsave("High clean net_clusters_Functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution


############################################ cl2_55_75

# Create the bar plot with dodge placement and fill by function
my_plot = ggplot(cl2_data, aes(x = Cluster, y = OTU_abundance, fill = Function)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Functional frequencies in medium cleanliness network", x = "Cluster", y = "Functional Frequency") +
  theme_classic()+
  #scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 1))
  scale_fill_manual(values = colors[match(cl2_data$Function, unique(cl2_data$Function))])+
  theme(legend.position = "right",           # Move legend to the right side
        legend.title = element_text(face = "bold"),  # Bold title
        legend.text = element_text(size = 7),legend.key.size = unit(0.3, "cm")) 
plot(my_plot)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Medium clean net_clusters_Functional frequencies.pdf", width = 11, height = 8.5, units = "in")  # Set maximum paper size (adjust units as needed)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  # Save as JPEG with maximum paper size (adjust units as needed)
  ggsave("Medium clean net_clusters_Functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution


############################################ cl3_35_55

# Create the bar plot with dodge placement and fill by function
my_plot = ggplot(cl3_data, aes(x = Cluster, y = OTU_abundance, fill = Function)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Functional frequencies in low cleanliness network", x = "Cluster", y = "Functional Frequency") +
  theme_classic() +
  #scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 1))
  scale_fill_manual(values = colors[match(cl3_data$Function, unique(cl3_data$Function))])+
  theme(legend.position = "right",           # Move legend to the right side
        legend.title = element_text(face = "bold"),  # Bold title
        legend.text = element_text(size = 7),legend.key.size = unit(0.3, "cm")) 

plot(my_plot)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Low clean net_clusters_Functional frequencies.pdf", width = 11, height = 8.5, units = "in")  # Set maximum paper size (adjust units as needed)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  # Save as JPEG with maximum paper size (adjust units as needed)
  ggsave("Low clean net_clusters_Functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution


############################################ large

# Create the bar plot with dodge placement and fill by function
my_plot = ggplot(large_data, aes(x = Cluster, y = OTU_abundance, fill = Function)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Functional frequencies in large farm size network", x = "Cluster", y = "Functional Frequency") +
  theme_classic() +
  #scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 1))
  scale_fill_manual(values = colors[match(large_data$Function, unique(large_data$Function))])+
  theme(legend.position = "right",           # Move legend to the right side
        legend.title = element_text(face = "bold"),  # Bold title
        legend.text = element_text(size = 7),legend.key.size = unit(0.3, "cm")) 

plot(my_plot)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Large net_clusters_Functional frequencies.pdf", width = 11, height = 8.5, units = "in")  # Set maximum paper size (adjust units as needed)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  # Save as JPEG with maximum paper size (adjust units as needed)
  ggsave("Large net_clusters_Functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution


############################################ med_small

# Create the bar plot with dodge placement and fill by function
my_plot = ggplot(med_small_data, aes(x = Cluster, y = OTU_abundance, fill = Function)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Functional frequencies in med_small farm size network", x = "Cluster", y = "Functional Frequency") +
  theme_classic() +
  #scale_y_continuous(breaks = seq(0, max(data$different_OTUs_for_the_function), by = 1))
  scale_fill_manual(values = colors[match(med_small_data$Function, unique(med_small_data$Function))])+
  theme(legend.position = "right",           # Move legend to the right side
        legend.title = element_text(face = "bold"),  # Bold title
        legend.text = element_text(size = 7),legend.key.size = unit(0.3, "cm")) 

plot(my_plot)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  ggsave("Med-Small net_clusters_Functional frequencies.pdf", width = 11, height = 8.5, units = "in")  # Set maximum paper size (adjust units as needed)

theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +  # Adjust margins for landscape
  # Save as JPEG with maximum paper size (adjust units as needed)
  ggsave("Med-Small net_clusters_Functional frequencies.jpeg", width = 11, height = 8.5, units = "in", dpi = 350)  # Set high resolution

