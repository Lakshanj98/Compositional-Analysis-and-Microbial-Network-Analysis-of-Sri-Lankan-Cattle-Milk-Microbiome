# load the packages
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(microbiome)
library(NetCoMi)
library(WGCNA)
library(stringr)
library(limma)
setwd("D:/Research/R codes")



#Adding .qza files to phyloseq object
physeq=qza_to_phyloseq(
  features = "table.qza",
  taxonomy="taxonomy.qza",
  tree="rooted-tree.qza"
)

# write the OTU table into an excel using for the imageGP Faprotax
setwd("D:/Research/R codes/Functional Annotations/ImageGP_FAPROTAX")
write.csv(physeq@otu_table, file="otu_table.csv")
write.csv(physeq@tax_table, file="taxonomy_table.csv")

write.table(physeq@otu_table, file="otu_table.tsv", sep="\t", row.names=TRUE)
write.table(physeq@tax_table, file="taxonomy_table.tsv", sep="\t", row.names=TRUE)

otu_dframe <- as.data.frame(physeq@otu_table)
tax_dframe <- as.data.frame(physeq@tax_table)
write_tsv(otu_dframe, file="otu_table.tsv")
write_tsv(tax_dframe, file="taxonomy_table.tsv")
setwd("D:/Research/R codes")

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

# plotting the heatmap
rank_names(complete_physeq)
top20 <- names(sort(taxa_sums(complete_physeq), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, complete_physeq)
plot_heatmap(complete_physeq)

plot_heatmap(ps.top20,taxa.label="Genus", sample.label="cleanliness", sample.order="cleanliness")
plot_heatmap(ps.top20,taxa.label="Genus", sample.label="FarmSize", sample.order="FarmSize")
head(sample_data(complete_physeq))

# create a prevalence table
#Lets generate a prevelance table (number of samples each taxa occurs in) for each taxa.
prevelancedf = apply(X = otu_table(complete_physeq),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(complete_physeq),
                          tax_table(complete_physeq))
prevelancedf[1:10,]
write.csv(prevelancedf,'impData/prevelancedf.csv')
# this table contains prevalence and total abundance.
# It contains prevalence of OTUs and total abundance.

# Now lets investigate low prevalence/abundance phylum and subset them out.
taxa_to_remove =plyr::ddply(prevelancedf, "Phylum", function(df1){data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
write.csv(taxa_to_remove,'impData/prevelancephylums.csv')
phyla2Filter = c("Synergistetes", "OP8", "TM6",
                 "SR1","GN02","Fibrobacteres",
                 "Nitrospirae","WPS-2","Planctomycetes")
                  
# Filter entries with unidentified Phylum.
# unknown <- subset_taxa(physeq = tax_table(complete_physeq), !is.na(Phylum) & !Phylum %in% c("", "NA"))
phyfiltered = subset_taxa(complete_physeq, !Phylum %in% phyla2Filter)

phyfiltered

###########  Alpha Diversity   ########################################################3
# saving alpha diversity notched boxplots for both variables
selectvs =c("cleanliness","FarmSize")
for (var1 in selectvs){
  #setwd("D:/Research/R codes/impData")
  
  #File name for notched boxplot
  trfname= paste( "AD", var1,"notched_boxplots.pdf", sep="_")
  jpegname = paste( "AD", var1,"notched_boxplots.jpeg", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
  
  p <- p + geom_boxplot( notch= TRUE, alpha=0.2) # notched boxplot
  
  print(plot(p))
  # Saving the pdf
  dev.off()
  
  # generate a high quality 350 ppi image
  jpeg(filename = jpegname,
       width = 2337, height = 1653, units = "px", pointsize = 12,
       quality = 75, bg = "white", res = 350)
  plot(p)
  dev.off()
  
  
  #File name for boxplot
  trfname= paste( "AD", var1,"boxplots.pdf", sep="_")
  jpegname = paste( "AD", var1,"boxplots.jpeg", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
  
  p <- p + geom_boxplot( alpha=0.2) 
  
  print(plot(p))
  # Saving the pdf
  dev.off()
  
  # generate a high quality 350 ppi image
  jpeg(filename = jpegname,
       width = 2337, height = 1653, units = "px", pointsize = 12,
       quality = 75, bg = "white", res = 350)
  plot(p)
  dev.off()
  
  #File name for violinplot
  trfname= paste( "AD", var1,"violinplots.pdf", sep="_")
  jpegname = paste( "AD", var1,"violinplots.jpeg", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
  
  p <- p + geom_violin(  alpha=0.2) # violinplot
  
  print(plot(p))
  # Saving the pdf
  dev.off()
  
  # generate a high quality 350 ppi image
  jpeg(filename = jpegname,
       width = 2337, height = 1653, units = "px", pointsize = 12,
       quality = 75, bg = "white", res = 350)
  plot(p)
  dev.off()
}

#generating richness values to a table
richtable =estimate_richness(phyfiltered)
# saving the table
write.csv(richtable,'ad_estimates.csv')


# list to combine wilcox test outputs
pwclist = list()

for (var1 in selectvs){
  pwc= pairwise.wilcox.test(richtable$Chao1, sample_data(phyfiltered)[[var1]],p.adj = "bonf")
  pwclist= c(pwclist,pwc)
}

for (var1 in selectvs){
  pwc= pairwise.wilcox.test(richtable$Observed, sample_data(phyfiltered)[[var1]],p.adj = "bonf")
  pwclist= c(pwclist,pwc)
}

for (var1 in selectvs){
  pwc= pairwise.wilcox.test(richtable$Shannon, sample_data(phyfiltered)[[var1]],p.adj = "bonf")
  pwclist= c(pwclist,pwc)
}

for (var1 in selectvs){
  pwc= pairwise.wilcox.test(richtable$Simpson, sample_data(phyfiltered)[[var1]],p.adj = "bonf")
  pwclist= c(pwclist,pwc)
  # capturing the combined outputs and saving in a text file
  chars= capture.output(print(pwclist))
  writeLines(chars, con = file("_filtered_bonfadjustedpvalue_four indexes.txt"))
  # closing the opened file connection
  close(file("_filtered_bonfadjustedpvalue_four indexes.txt"))
}

################# Beta diversity measures   #####################################

for (i in selectvs){
  #setwd("D:/Research/R codes/impData")
  #File name
  trfname= paste("Edited",i,"network.pdf", sep="_")
  jpegname= paste("Edited",i,"network.jpeg", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_network(igbray, phyfiltered, color=i, line_weight=0.4, label = i))
  # Saving the pdf
  dev.off()
  
  # generate a high quality 350 ppi image
  jpeg(filename = jpegname,
       width = 2337, height = 1653, units = "px", pointsize = 12,
       quality = 75, bg = "white", res = 350)
  #plotting
  print(plot_network(igbray, phyfiltered, color=i, line_weight=0.4, label = i))
  # Saving the image
  dev.off()
}

################################# Ordination  plots ######
# Calculate distances
DistBC = distance(phyfiltered, method = "bray")
DistUF = distance(phyfiltered, method = "wUniFrac")

# Calculating the ordinations
ordBC = ordinate(phyfiltered, method = "PCoA", distance = DistBC)
ordUF = ordinate(phyfiltered, method = "PCoA", distance = DistUF)

# Plotting scree plots
plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_scree(ordUF, "Scree Plot: Weighted UniFrac MDS")

# For bray-curtis
for (i in selectvs){
  #setwd("D:/Research/R codes/impData")
  #File name
  trfname= paste("Edited_Bray-curtis_ordinationplot",i,".pdf", sep="_")
  jpegname= paste("Edited_Bray-curtis_ordinationplot",i,".jpeg", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_ordination(phyfiltered, ordBC, color=i) +ggtitle("PCoA: Bray-Curtis"))
  # Saving the pdf
  dev.off()
  
  # generate a high quality 350 ppi image
  jpeg(filename = jpegname,
       width = 2337, height = 1653, units = "px", pointsize = 12,
       quality = 75, bg = "white", res = 350)
  #plotting
  print(plot_ordination(phyfiltered, ordBC, color=i) +ggtitle("PCoA: Bray-Curtis"))
  # Saving the pdf
  dev.off()
}



# For unifrac
for (i in selectvs){
  #setwd("D:/Research/R codes/impData")
  #File name
  trfname= paste("Edited_Unifrac_ordinationplot",i,".pdf", sep="_")
  jpegname= paste("Edited_Unifrac_ordinationplot",i,".jpeg", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_ordination(phyfiltered, ordUF, color=i, label=i) +ggtitle("PCoA: Weigthed Unifrac"))
  # Saving the pdf
  dev.off()
  
  # generate a high quality 350 ppi image
  jpeg(filename = jpegname,
       width = 2337, height = 1653, units = "px", pointsize = 12,
       quality = 75, bg = "white", res = 350)
  #plotting
  print(plot_ordination(phyfiltered, ordUF, color=i, label=i) +ggtitle("PCoA: Weigthed Unifrac"))
  # Saving the pdf
  dev.off()
}

################################      add the bar plots from line 423
