# load the packages
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(microbiome)
library(NetCoMi)
library(WGCNA)
library(stringr)
library(limma)
#install.packages("MatrixExtra")
library(MatrixExtra)
setwd("D:/Research/R codes")

physeq=qza_to_phyloseq(
  features = "table.qza",
  taxonomy="taxonomy.qza"
)
#Reading the meta data
metx=read.delim("metadata_table.txt",sep=",",header=T,row.names=sample_names(physeq))
head(metx)

#Meta data as sample data
metx=metx[,-1]
head(metx)
metx=sample_data(metx)
metx

#Merging the sample
complete_physeq=merge_phyloseq(physeq,metx)

complete_physeq

#complete_physeq
xa=sample_data(complete_physeq)
xb=tax_table(complete_physeq)
xc=otu_table(complete_physeq)

# checking whether there are incomplete samples
# and dropping if present.
ps_drop_incomplete <- function(ps, vars = NA, verbose = "max") {
  df <- sample_data(ps)
  if (identical(vars, NA)) vars <- phyloseq::sample_variables(ps)
  df_sub <- df[, vars, drop = FALSE]
  df_sub <- df_sub[stats::complete.cases(df_sub), , drop = FALSE]
  if (!isFALSE(verbose)) {
    incomplete <- nrow(df) - nrow(df_sub)
    if (incomplete > 0 || identical(verbose, "max")) {
      message("Dropping samples with missings: ", incomplete)
    }
  }
  if (identical(verbose, "max")) {
    for (v in vars) {
      n_missings <- sum(is.na(df[[v]]))
      if (n_missings > 0) message(v, " has NAs: ", n_missings)
    }
  }
  keepers <- rownames(df_sub)
  phyloseq::prune_samples(samples = keepers, x = ps)
}

complete_physeq=ps_drop_incomplete(complete_physeq)

# create a prevalence table  ----------------------------------------------
prevelancedf = apply(X = otu_table(complete_physeq),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(complete_physeq),
                          tax_table(complete_physeq))
prevelancedf[1:10,]
write.csv(prevelancedf,'prevelancedf.csv')
# Now lets investigate low prevelance/abundance phylum and subset them out.
taxa_to_remove =plyr::ddply(prevelancedf, "Phylum", function(df1){data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
write.csv(taxa_to_remove,'prevelancephylums.csv')

# filtering the taxa based on prevalance
phyla2Filter = c(#"Armatimonadetes" ,
                 "Synergistetes", "OP8", "TM6",
                 "SR1","GN02","Fibrobacteres",
                 "Nitrospirae","WPS-2","Planctomycetes")

# Filter entries with unidentified Phylum.
# unknown <- subset_taxa(physeq = tax_table(complete_physeq), !is.na(Phylum) & !Phylum %in% c("", "NA"))
phyfiltered = subset_taxa(complete_physeq, !Phylum %in% phyla2Filter)


phyfiltered



tax = tax_table(phyfiltered)
write.csv(tax, "taxonomy table.csv")

otu= otu_table(phyfiltered)
write.csv(otu, "otu table phy.csv")

write.table(otu, "otu table phy.tsv", sep="\t")

# transforming to obtain relative abundances
transformed = transform_sample_counts(phyfiltered, function(x) x / sum(x))

# # code part 
# 
# comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
# agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)
# 
# # filter out the genus column with NA values
# genus_filter = c("","NA")
# phyfiltered_g_filtered = subset_taxa(phyfiltered, !Genus %in% genus_filter)
# phyfiltered_g_filtered

#aggregating taxa in the genus level
phyfiltered_glom <- tax_glom(phyfiltered, 'Genus')

tax = tax_table(phyfiltered_glom)
write.csv(tax, "taxonomy table_glom.csv")
#complete_agg
tax_table(phyfiltered_glom)
final_all_data <- aggregate_taxa(phyfiltered_glom, 'Genus')
final_all_data

tax = tax_table(final_all_data)
write.csv(tax, "taxonomy table_agg.csv")

otu = otu_table(final_all_data)
write.csv(otu, "OTU table_agg.csv")

saveRDS(final_all_data, "networks/Sparcc/final_all_data.rds")

#saveRDS(phyfiltered, "networks/Sparcc/phyfiltered.rds")

# Co-occurrence network building ###############################################
################################################################################

### ---------------------- SparCC method ---------------------------------------- ###
#all_data <- readRDS("networks/Sparcc/phyfiltered.rds")
all_data <- readRDS("networks/Sparcc/final_all_data.rds")

#full_network <- prune_taxa(names(sort(taxa_sums(all_data),TRUE)[1:50]), all_data)

taxtab <- all_data@tax_table@.Data
phyla <- taxtab[, "Genus"]
plot_heatmap(full_network, taxa.label = phyla)
full_network


#SparCC
full_net_sparcc <- netConstruct(all_data,
                         verbose = 3,
                         filtTax = "highestFreq",
                         filtTaxPar = list(highestFreq = 100),
                         # filtTax = "none",
                         # filtTaxPar = "totalReads",
                         filtSamp = "totalReads",
                         filtSampPar = list(totalReads = 1000),
                         zeroMethod = "none", normMethod = "none",
                         measure = "sparcc",
                         sparsMethod = "threshold", thresh = 0.4,
                         #sparsMethod	
                         #the method used for sparsification (selected edges that are connected in the network). 
                         #    "threshold"
                         #Selected are taxa pairs with an absolute association/dissimilarity greater than or equal to the threshold defined via thresh.
                         dissFunc = "signed", 
                         seed = 123456)

saveRDS(full_net_sparcc, "networks/Sparcc/full_network.rds")


props_full_sparcc <- netAnalyze(full_net_sparcc, 
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector", 
                         hubQuant = 0.95)

saveRDS(props_full_sparcc, "networks/Sparcc/full_network_Sparcc analysis.rds")


# write the network summary to a csv file
write.csv(full_net_sparcc$edgelist1, file ="networks/Sparcc/full_net/full_net_sparcc.csv")
write.csv(full_net_plot_sparcc$labels$labels1, file ="networks/Sparcc/full_net/full_net_labels_sparcc.csv")
write.csv(props_full_sparcc$clustering$clust1, file ="networks/Sparcc/full_net/full_net_clusters_sparcc.csv")

#?summary.microNetProps, plot.microNetProps
net.summary <- summary(props_full_sparcc)
net.summary <- as.data.frame(net.summary)

# write the network summary to a csv file
write.csv(net.summary, file ="networks/Sparcc/full_network_sparcc.csv")

pdf("networks/Sparcc/full_network_sparcc.pdf", width = 30, height = 30)
full_net_plot_sparcc <- plot(props_full_sparcc,
                      shortenLabels = "none",
                      # labelLength = 16,
                      # charToRm = "g__",
                      labelScale = FALSE,
                      rmSingles = "all",
                      labels = phyla,
                      nodeSize = "eigenvector",
                      nodeColor = "cluster",
                      hubBorderCol = "blue",
                      cexNodes = 1,
                      cexLabels = 1.5,
                      labelFont = 2, 
                      cexHubLabels = 2, 
                      hubLabelFont = 2,
                      edgeWidth = 1,
                      highlightHubs = TRUE,
                      cexHubs = 1.5,
                      # cexHubLabels = 2,
                      title1 = "Full Network on Genus level with SparCC Method", 
                      showTitle = TRUE,
                      cexTitle = 2)

dev.off()

# generate a high quality 350 ppi image
jpeg(filename = "networks/Sparcc/full_network_sparcc.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)

full_net_plot_sparcc <- plot(props_full_sparcc,
                      shortenLabels = "none",
                      # labelLength = 16,
                      # charToRm = "g__",
                      labelScale = FALSE,
                      rmSingles = "all",
                      labels = phyla,
                      nodeSize = "eigenvector",
                      nodeColor = "cluster",
                      hubBorderCol = "blue",
                      cexNodes = 1,
                      cexLabels = 1.5,
                      labelFont = 2, 
                      cexHubLabels = 2, 
                      hubLabelFont = 2,
                      edgeWidth = 1,
                      highlightHubs = TRUE,
                      cexHubs = 1.5,
                      # cexHubLabels = 2,
                      title1 = "Full Network on Genus level with SparCC Method", 
                      showTitle = TRUE,
                      cexTitle = 2)
dev.off()

saveRDS(full_net_plot_sparcc, "networks/Sparcc/full_network_Sparcc plot.rds")

# regenerating the plot for extracting network details.
full_net_plot_sparcc <- plot(props_full_sparcc,
                             shortenLabels = "none",
                             # labelLength = 16,
                             # charToRm = "g__",
                             labelScale = FALSE,
                             rmSingles = "all",
                             #labels = phyla,
                             nodeSize = "eigenvector",
                             nodeColor = "cluster",
                             hubBorderCol = "blue",
                             cexNodes = 1,
                             cexLabels = 1.5,
                             labelFont = 2, 
                             cexHubLabels = 2, 
                             hubLabelFont = 2,
                             edgeWidth = 1,
                             highlightHubs = TRUE,
                             cexHubs = 1.5,
                             # cexHubLabels = 2,
                             title1 = "Full Network on Genus level with SparCC Method", 
                             showTitle = TRUE,
                             cexTitle = 2)

# extracting the network details

from <- full_net_plot_sparcc$labels$labels1[full_net_plot_sparcc$q1$Edgelist$from]
to <- full_net_plot_sparcc$labels$labels1[full_net_plot_sparcc$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (full_net_sparcc$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1
  }
  else if (full_net_sparcc$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1
  }
}

edges <- data.frame(row.names = c(1:length(full_net_plot_sparcc$q1$Edgelist$from)))
edges$from <- full_net_plot_sparcc$q1$Edgelist$from
edges$to <- full_net_plot_sparcc$q1$Edgelist$to
edges$weight <- full_net_plot_sparcc$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "networks/Sparcc/full_net/full_network_Sparcc edge data.csv", row.names = FALSE)

hubs <- props_full_sparcc$hubs$hubs1
write(hubs, "networks/Sparcc/full_net/full_network_Sparcc Hubs.txt")

node.lables <- full_net_plot_sparcc$labels$labels1
clust <- props_full_sparcc$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_full_sparcc$centralities$degree1[nodes$lable]
evs = props_full_sparcc$centralities$eigenv1[nodes$lable]
betweennesses = props_full_sparcc$centralities$between1[nodes$lable]
closenesses = props_full_sparcc$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses



write.csv(nodes, file = "networks/Sparcc/full_net/full_network_Sparcc node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "networks/Sparcc/full_net/full_network_Sparcc hub data.csv", row.names = FALSE)


### --------------------------------------------------------------------------
# Create phyloseq objects to all levels in variables.
cl1_75_95 = subset_samples(all_data, cleanliness == "cl1_75_95" )
cl2_55_75 = subset_samples(all_data, cleanliness == "cl2_55_75" )
cl3_35_55 = subset_samples(all_data, cleanliness == "cl3_35_55" )
Med_Small = subset_samples(all_data, FarmSize == "Med_Small" )
Large = subset_samples(all_data, FarmSize == "Large" )

# create network for cleanliness  == cl1_75_95  ####################################

#network_cl1_75_95 <- prune_taxa(names(sort(taxa_sums(cl1_75_95),TRUE)[1:50]), cl1_75_95)
#plot_heatmap(top_Early)
taxtab <- cl1_75_95@tax_table@.Data
phyla <- taxtab[, "Genus"]

# #aggregating taxa in the genus level
# phyfiltered_glom <- tax_glom(cl1_75_95, 'Genus')
# 
# #complete_agg
# cl1_75_95 <- aggregate_taxa(phyfiltered_glom, 'Genus')

#SparCC
net_cl1_75_95_sparcc <- netConstruct(cl1_75_95,
                              verbose = 3,
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              # filtTax = "none",
                              # filtTaxPar = "totalReads",
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              zeroMethod = "none", normMethod = "none",
                              measure = "sparcc",
                              sparsMethod = "threshold", thresh = 0.4,
                              dissFunc = "signed", 
                              seed = 123456)

saveRDS(net_cl1_75_95_sparcc, "networks/Sparcc/net_cl1_75_95_sparcc.rds")

props_cl1_75_95_sparcc <- netAnalyze(net_cl1_75_95_sparcc, 
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector", 
                              hubQuant = 0.95)

saveRDS(props_cl1_75_95_sparcc, "networks/Sparcc/net_cl1_75_95_sparcc_analysis.rds")

# write the network details to a csv
write.csv(net_cl1_75_95_sparcc$edgelist1, file ="networks/Sparcc/cl1_75_95/cl1_75_95_sparcc.csv")
write.csv(cl1_75_95_plot_sparcc$labels$labels1, file ="networks/Sparcc/cl1_75_95/cl1_75_95_labels_sparcc.csv")
write.csv(props_cl1_75_95_sparcc$clustering$clust1, file ="networks/Sparcc/cl1_75_95/cl1_75_95_clusters_sparcc.csv")

#?summary.microNetProps
net.summary <- summary(props_cl1_75_95_sparcc)
net.summary <- as.data.frame(net.summary)

pdf("networks/Sparcc/Cleanliness_cl1_75_95_sparcc_net.pdf", width = 30, height = 30)
cl1_75_95_plot_sparcc <- plot(props_cl1_75_95_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "Cleanliness High on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)

dev.off()

saveRDS(cl1_75_95_plot_sparcc, "networks/Sparcc/cl1_75_95_sparcc plot.rds")

jpeg(filename = "networks/Sparcc/Cleanliness_cl1_75_95_sparcc_net.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)

cl1_75_95_plot_sparcc <- plot(props_cl1_75_95_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "Cleanliness High on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)
dev.off()

# regenerating the plot for extracting network details.
cl1_75_95_plot_sparcc <- plot(props_cl1_75_95_sparcc,
                              shortenLabels = "none",
                              # labelLength = 16,
                              # charToRm = "g__",
                              labelScale = FALSE,
                              #labels = phyla,
                              rmSingles = "all",
                              nodeSize = "eigenvector",
                              nodeColor = "cluster",
                              hubBorderCol = "blue",
                              cexNodes = 1,
                              cexLabels = 1.5,
                              labelFont = 2, 
                              cexHubLabels = 2, 
                              hubLabelFont = 2,
                              edgeWidth = 1,
                              highlightHubs = TRUE,
                              cexHubs = 1.5,
                              # cexHubLabels = 2,
                              title1 = "Cleanliness High on Genus level with SparCC Method", 
                              showTitle = TRUE,
                              cexTitle = 2)

# extracting the network details
from <- cl1_75_95_plot_sparcc$labels$labels1[cl1_75_95_plot_sparcc$q1$Edgelist$from]
to <- cl1_75_95_plot_sparcc$labels$labels1[cl1_75_95_plot_sparcc$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_cl1_75_95_sparcc$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1
  }
  else if (net_cl1_75_95_sparcc$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1
  }
}

edges <- data.frame(row.names = c(1:length(cl1_75_95_plot_sparcc$q1$Edgelist$from)))
edges$from <- cl1_75_95_plot_sparcc$q1$Edgelist$from
edges$to <- cl1_75_95_plot_sparcc$q1$Edgelist$to
edges$weight <- cl1_75_95_plot_sparcc$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "networks/Sparcc/cl1_75_95/cl1_75_95_Sparcc edge data.csv", row.names = FALSE)

hubs <- props_cl1_75_95_sparcc$hubs$hubs1
write(hubs, "networks/Sparcc/cl1_75_95/cl1_75_95_Sparcc Hubs.txt")

node.lables <- cl1_75_95_plot_sparcc$labels$labels1
clust <- props_cl1_75_95_sparcc$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_cl1_75_95_sparcc$centralities$degree1[nodes$lable]
evs = props_cl1_75_95_sparcc$centralities$eigenv1[nodes$lable]
betweennesses = props_cl1_75_95_sparcc$centralities$between1[nodes$lable]
closenesses = props_cl1_75_95_sparcc$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses



write.csv(nodes, file = "networks/Sparcc/cl1_75_95/cl1_75_95_Sparcc node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "networks/Sparcc/cl1_75_95/cl1_75_95_Sparcc hub data.csv", row.names = FALSE)

# create network for cleanliness  == cl2_55_75 -----------------------------------
taxtab <- cl2_55_75@tax_table@.Data
phyla <- taxtab[, "Genus"]

#SparCC
net_cl2_55_75_sparcc <- netConstruct(cl2_55_75,
                              verbose = 3,
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              # filtTax = "none",
                              # filtTaxPar = "totalReads",
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              zeroMethod = "none", normMethod = "none",
                              measure = "sparcc",
                              sparsMethod = "threshold", thresh = 0.4,
                              dissFunc = "signed", 
                              seed = 123456)


saveRDS(net_cl2_55_75_sparcc, "networks/Sparcc/net_cl2_55_75_sparcc.rds")

props_cl2_55_75_sparcc <- netAnalyze(net_cl2_55_75_sparcc, 
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector", 
                              hubQuant = 0.95)

saveRDS(props_cl2_55_75_sparcc, "networks/Sparcc/net_cl2_55_75_sparcc_analysis.rds")

# write the network details to a csv
write.csv(net_cl2_55_75_sparcc$edgelist1, file ="networks/Sparcc/cl2_55_75/cl2_55_75_sparcc.csv")
write.csv(cl2_55_75_plot_sparcc$labels$labels1, file ="networks/Sparcc/cl2_55_75/cl2_55_75_labels_sparcc.csv")
write.csv(props_cl2_55_75_sparcc$clustering$clust1, file ="networks/Sparcc/cl2_55_75/cl2_55_75_clusters_sparcc.csv")


#?summary.microNetProps
net.summary <- summary(props_cl2_55_75_sparcc)
net.summary <- as.data.frame(net.summary)


pdf("networks/Sparcc/Cleanliness_cl2_55_75_sparcc_net.pdf", width = 30, height = 30)

cl2_55_75_plot_sparcc <- plot(props_cl2_55_75_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "Cleanliness Medium on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)

dev.off()

saveRDS(cl2_55_75_plot_sparcc, "networks/Sparcc/cl2_55_75_sparcc plot.rds")

jpeg(filename = "networks/Sparcc/Cleanliness_cl2_55_75_sparcc_net.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)

cl2_55_75_plot_sparcc <- plot(props_cl2_55_75_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "Cleanliness Medium on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)
dev.off()

cl2_55_75_plot_sparcc <- plot(props_cl2_55_75_sparcc,
                              shortenLabels = "none",
                              # labelLength = 16,
                              # charToRm = "g__",
                              labelScale = FALSE,
                              #labels = phyla,
                              rmSingles = "all",
                              nodeSize = "eigenvector",
                              nodeColor = "cluster",
                              hubBorderCol = "blue",
                              cexNodes = 1,
                              cexLabels = 1.5,
                              labelFont = 2, 
                              cexHubLabels = 2, 
                              hubLabelFont = 2,
                              edgeWidth = 1,
                              highlightHubs = TRUE,
                              cexHubs = 1.5,
                              # cexHubLabels = 2,
                              title1 = "Cleanliness Medium on Genus level with SparCC Method", 
                              showTitle = TRUE,
                              cexTitle = 2)

# extracting the network details
from <- cl2_55_75_plot_sparcc$labels$labels1[cl2_55_75_plot_sparcc$q1$Edgelist$from]
to <- cl2_55_75_plot_sparcc$labels$labels1[cl2_55_75_plot_sparcc$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_cl2_55_75_sparcc$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1
  }
  else if (net_cl2_55_75_sparcc$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1
  }
}

edges <- data.frame(row.names = c(1:length(cl2_55_75_plot_sparcc$q1$Edgelist$from)))
edges$from <- cl2_55_75_plot_sparcc$q1$Edgelist$from
edges$to <- cl2_55_75_plot_sparcc$q1$Edgelist$to
edges$weight <- cl2_55_75_plot_sparcc$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "networks/Sparcc/cl2_55_75/cl2_55_75_Sparcc edge data.csv", row.names = FALSE)

hubs <- props_cl2_55_75_sparcc$hubs$hubs1
write(hubs, "networks/Sparcc/cl2_55_75/cl1_75_95_Sparcc Hubs.txt")

node.lables <- cl2_55_75_plot_sparcc$labels$labels1
clust <- props_cl2_55_75_sparcc$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_cl2_55_75_sparcc$centralities$degree1[nodes$lable]
evs = props_cl2_55_75_sparcc$centralities$eigenv1[nodes$lable]
betweennesses = props_cl2_55_75_sparcc$centralities$between1[nodes$lable]
closenesses = props_cl2_55_75_sparcc$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses



write.csv(nodes, file = "networks/Sparcc/cl2_55_75/cl2_55_75_Sparcc node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "networks/Sparcc/cl2_55_75/cl2_55_75_Sparcc hub data.csv", row.names = FALSE)


# create network for cleanliness  == cl3_35_55
#network_cl3_35_55 <- prune_taxa(names(sort(taxa_sums(cl3_35_55),TRUE)[1:50]), cl3_35_55)
taxtab <- cl3_35_55@tax_table@.Data
phyla <- taxtab[, "Genus"]

#SparCC
net_cl3_35_55_sparcc <- netConstruct(cl3_35_55,
                              verbose = 3,
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              # filtTax = "none",
                              # filtTaxPar = "totalReads",
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              zeroMethod = "none", normMethod = "none",
                              measure = "sparcc",
                              sparsMethod = "threshold", thresh = 0.4,
                              dissFunc = "signed", 
                              seed = 123456)

saveRDS(net_cl3_35_55_sparcc, "networks/Sparcc/net_cl3_35_55_sparcc.rds")

props_cl3_35_55_sparcc <- netAnalyze(net_cl3_35_55_sparcc, 
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector", 
                              hubQuant = 0.95)

saveRDS(props_cl3_35_55_sparcc, "networks/Sparcc/net_cl3_35_55_sparcc_analysis.rds")

# write the network details to a csv
write.csv(net_cl3_35_55_sparcc$edgelist1, file ="networks/Sparcc/cl3_35_55/cl3_35_55_sparcc.csv")
write.csv(cl3_35_55_plot_sparcc$labels$labels1, file ="networks/Sparcc/cl3_35_55/cl3_35_55_labels_sparcc.csv")
write.csv(props_cl3_35_55_sparcc$clustering$clust1, file ="networks/Sparcc/cl3_35_55/cl3_35_55_clusters_sparcc.csv")


#?summary.microNetProps
net.summary <- summary(props_cl3_35_55_sparcc)
net.summary <- as.data.frame(net.summary)

pdf("networks/Sparcc/Cleanliness_cl3_35_55_sparcc_net.pdf", width = 30, height = 30)

cl3_35_55_plot_sparcc <- plot(props_cl3_35_55_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "Cleanliness Low on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)

dev.off()

saveRDS(cl3_35_55_plot_sparcc, "networks/Sparcc/cl3_35_55_sparcc plot.rds")

jpeg(filename = "networks/Sparcc/Cleanliness_cl3_35_55_sparcc_net.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)

cl3_35_55_plot_sparcc <- plot(props_cl3_35_55_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "Cleanliness Low on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)
dev.off()

cl3_35_55_plot_sparcc <- plot(props_cl3_35_55_sparcc,
                              shortenLabels = "none",
                              # labelLength = 16,
                              # charToRm = "g__",
                              labelScale = FALSE,
                              #labels = phyla,
                              rmSingles = "all",
                              nodeSize = "eigenvector",
                              nodeColor = "cluster",
                              hubBorderCol = "blue",
                              cexNodes = 1,
                              cexLabels = 1.5,
                              labelFont = 2, 
                              cexHubLabels = 2, 
                              hubLabelFont = 2,
                              edgeWidth = 1,
                              highlightHubs = TRUE,
                              cexHubs = 1.5,
                              # cexHubLabels = 2,
                              title1 = "Cleanliness Low on Genus level with SparCC Method", 
                              showTitle = TRUE,
                              cexTitle = 2)

# extracting the network details
from <- cl3_35_55_plot_sparcc$labels$labels1[cl3_35_55_plot_sparcc$q1$Edgelist$from]
to <- cl3_35_55_plot_sparcc$labels$labels1[cl3_35_55_plot_sparcc$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_cl3_35_55_sparcc$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1
  }
  else if (net_cl3_35_55_sparcc$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1
  }
}

edges <- data.frame(row.names = c(1:length(cl3_35_55_plot_sparcc$q1$Edgelist$from)))
edges$from <- cl3_35_55_plot_sparcc$q1$Edgelist$from
edges$to <- cl3_35_55_plot_sparcc$q1$Edgelist$to
edges$weight <- cl3_35_55_plot_sparcc$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "networks/Sparcc/cl3_35_55/cl3_35_55_Sparcc edge data.csv", row.names = FALSE)

hubs <- props_cl3_35_55_sparcc$hubs$hubs1
write(hubs, "networks/Sparcc/cl3_35_55/cl3_35_55_Sparcc Hubs.txt")

node.lables <- cl3_35_55_plot_sparcc$labels$labels1
clust <- props_cl3_35_55_sparcc$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_cl3_35_55_sparcc$centralities$degree1[nodes$lable]
evs = props_cl3_35_55_sparcc$centralities$eigenv1[nodes$lable]
betweennesses = props_cl3_35_55_sparcc$centralities$between1[nodes$lable]
closenesses = props_cl3_35_55_sparcc$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses



write.csv(nodes, file = "networks/Sparcc/cl3_35_55/cl3_35_55_Sparcc node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "networks/Sparcc/cl3_35_55/cl3_35_55_Sparcc hub data.csv", row.names = FALSE)

# create network for FarmSize = Med_Small ------------------------------------
#network_Med_Small <- prune_taxa(names(sort(taxa_sums(Med_Small),TRUE)[1:50]), Med_Small)
#plot_heatmap(top_Early)

taxtab <- Med_Small@tax_table@.Data
phyla <- taxtab[, "Genus"]

#SparCC
net_Med_Small_sparcc <- netConstruct(Med_Small,
                              verbose = 3,
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              # filtTax = "none",
                              # filtTaxPar = "totalReads",
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              zeroMethod = "none", normMethod = "none",
                              measure = "sparcc",
                              sparsMethod = "threshold", thresh = 0.4,
                              dissFunc = "signed", 
                              seed = 123456)

saveRDS(net_Med_Small_sparcc, "networks/Sparcc/net_Med_Small_sparcc.rds")

props_Med_Small_sparcc <- netAnalyze(net_Med_Small_sparcc, 
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector", 
                              hubQuant = 0.95)

saveRDS(props_Med_Small_sparcc, "networks/Sparcc/net_Med_Small_sparcc_analysis.rds")

# write the network summary to a csv file
write.csv(net_Med_Small_sparcc$edgelist1, file ="networks/Sparcc/Med_small/Med_Small_sparcc.csv")
write.csv(Med_Small_plot_sparcc$labels$labels1, file ="networks/Sparcc/Med_small/Med_Small_labels_sparcc.csv")
write.csv(props_Med_Small_sparcc$clustering$clust1, file ="networks/Sparcc/Med_small/Med_Small_clusters_sparcc.csv")


#?summary.microNetProps
net.summary <- summary(props_Med_Small_sparcc)
net.summary <- as.data.frame(net.summary)

pdf("networks/Sparcc/FarmSize_Med_Small_sparcc_net.pdf", width = 30, height = 30)

Med_Small_plot_sparcc <- plot(props_Med_Small_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "FarmSize Med_Small on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)


dev.off()

saveRDS(Med_Small_plot_sparcc, "networks/Sparcc/Med_Small_sparcc plot.rds")

jpeg(filename = "networks/Sparcc/FarmSize_Med_Small_sparcc_net.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)

Med_Small_plot_sparcc <- plot(props_Med_Small_sparcc,
                       shortenLabels = "none",
                       # labelLength = 16,
                       # charToRm = "g__",
                       labelScale = FALSE,
                       labels = phyla,
                       rmSingles = "all",
                       nodeSize = "eigenvector",
                       nodeColor = "cluster",
                       hubBorderCol = "blue",
                       cexNodes = 1,
                       cexLabels = 1.5,
                       labelFont = 2, 
                       cexHubLabels = 2, 
                       hubLabelFont = 2,
                       edgeWidth = 1,
                       highlightHubs = TRUE,
                       cexHubs = 1.5,
                       # cexHubLabels = 2,
                       title1 = "FarmSize Med_Small on Genus level with SparCC Method", 
                       showTitle = TRUE,
                       cexTitle = 2)
dev.off()


Med_Small_plot_sparcc <- plot(props_Med_Small_sparcc,
                              shortenLabels = "none",
                              # labelLength = 16,
                              # charToRm = "g__",
                              labelScale = FALSE,
                              #labels = phyla,
                              rmSingles = "all",
                              nodeSize = "eigenvector",
                              nodeColor = "cluster",
                              hubBorderCol = "blue",
                              cexNodes = 1,
                              cexLabels = 1.5,
                              labelFont = 2, 
                              cexHubLabels = 2, 
                              hubLabelFont = 2,
                              edgeWidth = 1,
                              highlightHubs = TRUE,
                              cexHubs = 1.5,
                              # cexHubLabels = 2,
                              title1 = "FarmSize Med_Small on Genus level with SparCC Method", 
                              showTitle = TRUE,
                              cexTitle = 2)

# extracting the network details
from <- Med_Small_plot_sparcc$labels$labels1[Med_Small_plot_sparcc$q1$Edgelist$from]
to <- Med_Small_plot_sparcc$labels$labels1[Med_Small_plot_sparcc$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_Med_Small_sparcc$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1
  }
  else if (net_Med_Small_sparcc$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1
  }
}

edges <- data.frame(row.names = c(1:length(Med_Small_plot_sparcc$q1$Edgelist$from)))
edges$from <- Med_Small_plot_sparcc$q1$Edgelist$from
edges$to <- Med_Small_plot_sparcc$q1$Edgelist$to
edges$weight <- Med_Small_plot_sparcc$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "networks/Sparcc/Med_Small/Med_Small_Sparcc edge data.csv", row.names = FALSE)

hubs <- props_Med_Small_sparcc$hubs$hubs1
write(hubs, "networks/Sparcc/Med_Small/Med_Small_Sparcc Hubs.txt")

node.lables <- Med_Small_plot_sparcc$labels$labels1
clust <- props_Med_Small_sparcc$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_Med_Small_sparcc$centralities$degree1[nodes$lable]
evs = props_Med_Small_sparcc$centralities$eigenv1[nodes$lable]
betweennesses = props_Med_Small_sparcc$centralities$between1[nodes$lable]
closenesses = props_Med_Small_sparcc$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses



write.csv(nodes, file = "networks/Sparcc/Med_Small/Med_Small_Sparcc node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "networks/Sparcc/Med_Small/Med_Small_Sparcc hub data.csv", row.names = FALSE)

# create network for FarmSize = Large ################################################
#network_Large <- prune_taxa(names(sort(taxa_sums(Large),TRUE)[1:50]), Large)

taxtab <- Large@tax_table@.Data
phyla <- taxtab[, "Genus"]

#SparCC
net_Large_sparcc <- netConstruct(Large,
                          verbose = 3,
                          filtTax = "highestFreq",
                          filtTaxPar = list(highestFreq = 100),
                          # filtTax = "none",
                          # filtTaxPar = "totalReads",
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "sparcc",
                          sparsMethod = "threshold", thresh = 0.4,
                          dissFunc = "signed", 
                          seed = 123456)

saveRDS(net_Large_sparcc, "networks/Sparcc/net_Large_sparcc.rds")

props_Large_sparcc <- netAnalyze(net_Large_sparcc, 
                          clustMethod = "cluster_fast_greedy",
                          hubPar = "eigenvector", 
                          hubQuant = 0.95)

saveRDS(props_Large_sparcc, "networks/Sparcc/net_Large_sparcc_analysis.rds")

# write the network summary to a csv file
write.csv(net_Large_sparcc$edgelist1, file ="networks/Sparcc/Large/Large_sparcc.csv")
write.csv(Large_plot_sparcc$labels$labels1, file ="networks/Sparcc/Large/Large_labels_sparcc.csv")
write.csv(props_Large_sparcc$clustering$clust1, file ="networks/Sparcc/Large/Large_clusters_sparcc.csv")

#?summary.microNetProps
net.summary <- summary(props_Large_sparcc)
net.summary <- as.data.frame(net.summary)

pdf("networks/Sparcc/FarmSize_Large_sparcc_net.pdf", width = 30, height = 30)

Large_plot_sparcc <- plot(props_Large_sparcc,
                   shortenLabels = "none",
                   # labelLength = 16,
                   # charToRm = "g__",
                   labelScale = FALSE,
                   labels = phyla,
                   rmSingles = "all",
                   nodeSize = "eigenvector",
                   nodeColor = "cluster",
                   hubBorderCol = "blue",
                   cexNodes = 1,
                   cexLabels = 1.5,
                   labelFont = 2, 
                   cexHubLabels = 2, 
                   hubLabelFont = 2,
                   edgeWidth = 1,
                   highlightHubs = TRUE,
                   cexHubs = 1.5,
                   # cexHubLabels = 2,
                   title1 = "FarmSize Large on Genus level with SparCC Method", 
                   showTitle = TRUE,
                   cexTitle = 2)

dev.off()

saveRDS(Large_plot_sparcc, "networks/Sparcc/Large_sparcc plot.rds")

jpeg(filename = "networks/Sparcc/FarmSize_Large_sparcc_net.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)

Large_plot_sparcc <- plot(props_Large_sparcc,
                   shortenLabels = "none",
                   # labelLength = 16,
                   # charToRm = "g__",
                   labelScale = FALSE,
                   labels = phyla,
                   rmSingles = "all",
                   nodeSize = "eigenvector",
                   nodeColor = "cluster",
                   hubBorderCol = "blue",
                   cexNodes = 1,
                   cexLabels = 1.5,
                   labelFont = 2, 
                   cexHubLabels = 2, 
                   hubLabelFont = 2,
                   edgeWidth = 1,
                   highlightHubs = TRUE,
                   cexHubs = 1.5,
                   # cexHubLabels = 2,
                   title1 = "FarmSize Large on Genus level with SparCC Method", 
                   showTitle = TRUE,
                   cexTitle = 2)
dev.off()


Large_plot_sparcc <- plot(props_Large_sparcc,
                              shortenLabels = "none",
                              # labelLength = 16,
                              # charToRm = "g__",
                              labelScale = FALSE,
                              #labels = phyla,
                              rmSingles = "all",
                              nodeSize = "eigenvector",
                              nodeColor = "cluster",
                              hubBorderCol = "blue",
                              cexNodes = 1,
                              cexLabels = 1.5,
                              labelFont = 2, 
                              cexHubLabels = 2, 
                              hubLabelFont = 2,
                              edgeWidth = 1,
                              highlightHubs = TRUE,
                              cexHubs = 1.5,
                              # cexHubLabels = 2,
                              title1 = "FarmSize Large on Genus level with SparCC Method", 
                              showTitle = TRUE,
                              cexTitle = 2)

# extracting the network details
from <- Large_plot_sparcc$labels$labels1[Large_plot_sparcc$q1$Edgelist$from]
to <- Large_plot_sparcc$labels$labels1[Large_plot_sparcc$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_Large_sparcc$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1
  }
  else if (net_Large_sparcc$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1
  }
}

edges <- data.frame(row.names = c(1:length(Large_plot_sparcc$q1$Edgelist$from)))
edges$from <- Large_plot_sparcc$q1$Edgelist$from
edges$to <- Large_plot_sparcc$q1$Edgelist$to
edges$weight <- Large_plot_sparcc$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "networks/Sparcc/Large/Large_Sparcc edge data.csv", row.names = FALSE)

hubs <- props_Large_sparcc$hubs$hubs1
write(hubs, "networks/Sparcc/Large/Large_Sparcc Hubs.txt")

node.lables <- Large_plot_sparcc$labels$labels1
clust <- props_Large_sparcc$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_Large_sparcc$centralities$degree1[nodes$lable]
evs = props_Large_sparcc$centralities$eigenv1[nodes$lable]
betweennesses = props_Large_sparcc$centralities$between1[nodes$lable]
closenesses = props_Large_sparcc$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses



write.csv(nodes, file = "networks/Sparcc/Large/Large_Sparcc node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "networks/Sparcc/Large/Large_Sparcc hub data.csv", row.names = FALSE)



           

