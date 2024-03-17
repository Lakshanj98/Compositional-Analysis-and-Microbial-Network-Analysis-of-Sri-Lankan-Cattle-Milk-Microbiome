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

################################    Relative abundances in milk microbiota

# generating the basic plot
top100 <- names(sort(taxa_sums(phyfiltered), decreasing=TRUE))[1:100]
# prune_taxa: An S4 Generic method for removing (pruning) unwanted OTUs/taxa 
# from phylogenetic objects, including phylo-class trees, as well as native phyloseq package objects. 
ps.top100 <- prune_taxa(top100, phyfiltered)
plot_bar(ps.top100, fill="Order")

# generate a high quality 350 ppi image
jpeg(filename = "Top100_orders with abundance.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot_bar(ps.top100, fill="Order")
# Saving the pdf
dev.off()

# converting to RA values
transtot <- transform_sample_counts(phyfiltered, function(x) x / sum(x))
ps.top100 <- prune_taxa(top100, transtot)
plot_bar(ps.top100, fill="Order")+ylab("Relative abundance")
# generate a high quality 350 ppi image
jpeg(filename = "Top100_orders with relative abundance.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot_bar(ps.top100, fill="Order")+ylab("Relative abundance")
# Saving the pdf
dev.off()



# defining a taxa level
t_level = "Phylum"

# converting to RA values
transtot <- transform_sample_counts(phyfiltered, function(x) x / sum(x))
#aggregating taxa (aggregating RA at genus level )
totGlommed = tax_glom(transtot, t_level, NArm=TRUE)

top20 <- names(sort(taxa_sums(totGlommed), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, totGlommed)

# filtering based on the percentage
cbrfr = filter_taxa(totGlommed, function(x)  max(x) > 0.01, TRUE)
# use the following plot margin for fat bars
plot_bar(cbrfr, fill=t_level)+theme(  plot.margin = margin(2, 9, 2, 9, "cm"))+ ylab("Relative abundance")
plot_bar(cbrfr, fill=t_level)+ylab("Relative abundance")

# generate a high quality 350 ppi image
jpeg(filename = "Phyla greater than 0.01 with relative abundance.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot_bar(cbrfr, fill=t_level)+ylab("Relative abundance")
# Saving the pdf
dev.off()


plot_bar(ps.top20, fill=t_level)+ylab("Relative abundance")

# generate a high quality 350 ppi image
jpeg(filename = "Top20_Phylumss with relative abundance.jpeg",
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot_bar(ps.top20, fill=t_level)+ylab("Relative abundance")
# Saving the pdf
dev.off()



totmelted <- psmelt(totGlommed)

# concatenating string for csv file name
#     setwd("D:/Research/R codes/impData")
# writing relative abundance values of generas
csvfname= paste(t_level, "_tot_RA.csv", sep="_")
write.csv(totmelted, csvfname)


setwd("D:/Research/R codes")

# Generating phylogenetic trees

# not sure section
# This method merges species that have the same taxonomy at a certain taxonomic rank.
phylumGlommed = tax_glom(transformed, t_level, NArm=FALSE)

# getting the top abundant taxa
topphy = names(sort(taxa_sums(phylumGlommed), TRUE)[1:30])
# prune_taxa: An S4 Generic method for removing (pruning) unwanted OTUs/taxa 
# from phylogenetic objects, including phylo-class trees, as well as native phyloseq package objects. 
toptax=prune_taxa(topphy,phylumGlommed)


# aggregating taxa
aggrphy <- aggregate_taxa(transformed, 'Phylum')
melted <- psmelt(toptax)
# Saving the otu table of the aggregated phyloseq object
#write_phyloseq(aggrphy, 'OTU', path = getwd())
#setwd("D:/Research/R codes/impData")
write.csv(melted,'GEnus_percentage.csv')

##  -------- generate trees with top 30 ---------------------------
# defining a taxa level
t_levellist = c("Species","Order","Genus", "Phylum")

# The list of variables
selectvs =c("cleanliness","FarmSize")

# cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.01, TRUE)
plot_tree(toptax, color="cleanliness", label.tips= "taxa_names", ladderize = TRUE, justify = "huha",nodelabf = nodeplotboot())
plot_tree(toptax, color="FarmSize", label.tips= "taxa_names", ladderize = TRUE, justify = "huha",nodelabf = nodeplotboot())


for(t_level in t_levellist){
  
  # Another way to aggregate data
  phylumGlommed = tax_glom(transformed, t_level, NArm=TRUE)
  
  # getting the top abundant taxa
  topphy = names(sort(taxa_sums(phylumGlommed), TRUE)[1:30])
  toptax=prune_taxa(topphy,phylumGlommed)
  
  # Plotting the phylogenetic trees for all the variables using a loop
  for (i in selectvs){
    #File name
    trfname= paste( "Top30_justified",t_level,"tree", i,".pdf", sep="_")
    jpegname= paste( "Top30_justified",t_level,"tree", i,".jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    #plotting
    print(plot_tree(toptax, color=i, label.tips= t_level,ladderize = "left", justify = "yes"))
    # Saving the pdf
    dev.off()
    
    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2337, height = 1653, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    print(plot_tree(toptax, color=i, label.tips= t_level,ladderize = "left", justify = "yes"))
    # Saving the pdf
    dev.off()
  }
  
  
}



# Saving the melted data file
phylummelted <- psmelt(phylumGlommed)

write.csv(phylummelted,'Order_glom_percentage.csv')


### Phylogenetic trees for different variable aggregates ###############
# different taxa ranks to iterate
taxal = c("Genus")

# defining a taxa level
t_levellist = c("Species","Order","Genus", "Phylum")

# The list of variables
selectvs =c("cleanliness","FarmSize")

# Looping through the selected sample variables
for (i in selectvs){
  # merge based on variables
  mergedphy =merge_samples(phyfiltered,group = i)
  
  # transform to relative abundance
  cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
  # Aggregation based on given taxa level
  for(tl in t_levellist){
    phylumGlommed = tax_glom(cbreedphy, tl, NArm=FALSE)
    
    # filtering taxa based on percentage
    #cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.02, TRUE)
    
    topphy = names(sort(taxa_sums(phylumGlommed), TRUE)[1:30])
    toptax=prune_taxa(topphy,phylumGlommed)
    
    #File name
    trfname= paste( "top30",tl,"tree", i,".pdf", sep="_")
    jpegname= paste( "top30",tl,"tree", i,".jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    #plotting
    print(plot_tree(toptax, color="Sample", size="abundance",label.tips= tl,base.spacing=0.04, ladderize = TRUE))
    
    # Saving the pdf
    dev.off()
    
    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2337, height = 1653, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    print(plot_tree(toptax, color=i, label.tips= t_level,ladderize = "left", justify = "yes"))
    # Saving the pdf
    dev.off()
  }
}

#################################       Bar Plots                ###############
#Variable aggregation
  # different taxa ranks to iterate
  selectvs =c("cleanliness","FarmSize")
  taxal = c("Species")             # for Species level
  
  
  # Looping through the selected sample variables
  for (i in selectvs){
    # merge based on variables
    mergedphy =merge_samples(phyfiltered,group = i)
    
    # transform to relative abundance
    cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
    # Aggregation based on given taxa level
    for(tl in taxal){
      phylumGlommed = tax_glom(cbreedphy, tl, NArm=FALSE)
      
      # merging the OTU table with TAXA names and saving as a CSV
      # converting each table to dataframes
      ttdf=data.frame(as(tax_table(phylumGlommed), "matrix"))
      # have to transpose the OTU table as merging automatically transposes the matrix
      otdf=as.data.frame(t(otu_table(phylumGlommed)))
      # merging
      mtable= merge.data.frame(ttdf,otdf, by=0,all.x = F)
      # saving as a CSV
      csname= paste( tl, i,"relative_abandance.csv", sep="_")
      write.csv(mtable,csname)
      
      # filtering taxa based on percentage
      cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.01, TRUE)
      
      #Plotting the stacked barplots
      #File name
      trfname= paste( tl, i,"barplot.pdf", sep="_")
      jpegname= paste( tl, i,"barplot.jpeg", sep="_")
      # Command to save as a pdf
      pdf(file=trfname, paper="a4r",width = 0, height = 0)
      # Plotting the stacked barplots
      p = plot_bar(subset_taxa(cbrfr,!is.na(Species)),fill = tl) + ylab("Relative abundance") + xlab(i)
      p = p + geom_bar(color="black",stat="identity", position="stack")
      # p=p+scale_fill_manual(values=c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
      #                                "#9467bd","#8c564b","#e377c2","#7f7f7f",
      #                                "#bcbd22","#17becf","#aec7e8","#ffbb78",
      #                                "#98df8a","#ff9896","#c5b0d5","#c49c94"))
      plot(p)
      # Saving the pdf
      dev.off()
      
      # generate a high quality 350 ppi image
      jpeg(filename = jpegname,
           width = 2337, height = 1653, units = "px", pointsize = 12,
           quality = 75, bg = "white", res = 350)
      #plotting
      plot(p)
      # Saving the pdf
      dev.off()
      
      
      # Heatmaps
      #File name
      trfname= paste( tl, i,"heatmap.pdf", sep="_")
      jpegname= paste( tl, i,"heatmap.jpeg", sep="_")
      # Command to save as a pdf
      pdf(file=trfname, paper="a4r",width = 0, height = 0)
      hm=plot_heatmap(cbrfr,method=NULL, taxa.label = tl,taxa.order=tl) + ylab("Relative abundance") + xlab(i) + scale_fill_gradient(name ="Relative abundance")
      plot(hm)
      # Saving the pdf
      dev.off()
      
      # generate a high quality 350 ppi image
      jpeg(filename = jpegname,
           width = 2337, height = 1653, units = "px", pointsize = 12,
           quality = 75, bg = "white", res = 350)
      #plotting
      plot(hm)
      # Saving the pdf
      dev.off()
    }
  }
  
################      For Genus Level

taxal = c("Genus")
# Looping through the selected sample variables
for (i in selectvs){
  # merge based on variables
  mergedphy =merge_samples(phyfiltered,group = i)
  
  # transform to relative abundance
  cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
  # Aggregation based on given taxa level
  for(tl in taxal){
    phylumGlommed = tax_glom(cbreedphy, tl, NArm=FALSE)
    
    # merging the OTU table with TAXA names and saving as a CSV
    # converting each table to dataframes
    ttdf=data.frame(as(tax_table(phylumGlommed), "matrix"))
    # have to transpose the OTU table as merging automatically transposes the matrix
    otdf=as.data.frame(t(otu_table(phylumGlommed)))
    # merging
    mtable= merge.data.frame(ttdf,otdf, by=0,all.x = F)
    # saving as a CSV
    csname= paste( tl, i,"relative_abandance.csv", sep="_")
    write.csv(mtable,csname)
    
    # filtering taxa based on percentage
    cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.01, TRUE)
    
    #Plotting the stacked barplots
    #File name
    trfname= paste( tl, i,"barplot.pdf", sep="_")
    jpegname= paste( tl, i,"barplot.jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    # Plotting the stacked barplots
    p = plot_bar(subset_taxa(cbrfr,!is.na(Genus)),fill = tl) + ylab("Relative abundance")+ xlab(i)
    p = p + geom_bar(color="black",stat="identity", position="stack")
    p=p+scale_fill_manual(values=c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
                                   "#9467bd","#8c564b","#e377c2","#7f7f7f",
                                   "#bcbd22","#17becf","#aec7e8","#ffbb78",
                                   "#98df8a","#ff9896","#c5b0d5","#c49c94",
                                   "#f7b6d2","#c7c7c7","#dbdb8d","#9edae5"))
    plot(p)
    # Saving the pdf
    dev.off()
    
    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2400, height = 2060, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    plot(p)
    # Saving the pdf
    dev.off()
    
    # plot in Rstudio to check graph
    plot(p)
    
    # Heatmaps
    #File name
    trfname= paste( tl, i,"heatmap.pdf", sep="_")
    jpegname= paste( tl, i,"heatmap.jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    hm=plot_heatmap(cbrfr,method=NULL, taxa.label = tl,taxa.order=tl)+ ylab("Relative abundance") + xlab(i)+ scale_fill_gradient(name ="Relative abundance")
    plot(hm)
    # Saving the pdf
    dev.off()

    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2337, height = 1653, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    plot(hm)
    # Saving the pdf
    dev.off()
  }
  
  
}




###############################   Order Level

#setwd("D:/Research/R codes")
#Variable aggregation
# different taxa ranks to iterate
taxal = c("Order")

# Looping through the selected sample variables
for (i in selectvs){
  # merge based on variables
  mergedphy =merge_samples(phyfiltered,group = i)
  
  # transform to relative abundance
  cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
  # Aggregation based on given taxa level
  for(tl in taxal){
    phylumGlommed = tax_glom(cbreedphy, tl, NArm=FALSE)
    
    # merging the OTU table with TAXA names and saving as a CSV
    # converting each table to dataframes
    ttdf=data.frame(as(tax_table(phylumGlommed), "matrix"))
    # have to transpose the OTU table as merging automatically transposes the matrix
    otdf=as.data.frame(t(otu_table(phylumGlommed)))
    # merging
    mtable= merge.data.frame(ttdf,otdf, by=0,all.x = F)
    # saving as a CSV
    csname= paste( tl, i,"relative_abandance.csv", sep="_")
    write.csv(mtable,csname)
    
    # filtering taxa based on percentage
    cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.01, TRUE)
    
    #Plotting the stacked barplots
    #File name
    trfname= paste( tl, i,"barplot.pdf", sep="_")
    jpegname= paste( tl, i,"barplot.jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    # Plotting the stacked barplots
    p = plot_bar(subset_taxa(cbrfr,!is.na(Order)),fill = tl) + ylab("Relative abundance")+ xlab(i)
    p = p + geom_bar(color="black",stat="identity", position="stack")
    p=p+scale_fill_manual(values=c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
                                   "#9467bd","#8c564b","#e377c2","#7f7f7f",
                                   "#bcbd22","#17becf","#aec7e8","#ffbb78",
                                   "#98df8a","#ff9896","#c5b0d5","#c49c94"))
    plot(p)
    # Saving the pdf
    dev.off()
    
    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2337, height = 1653, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    plot(p)
    # Saving the pdf
    dev.off()
    
    
    
    # Heatmaps
    #File name
    trfname= paste( tl, i,"heatmap.pdf", sep="_")
    jpegname= paste( tl, i,"heatmap.jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    hm=plot_heatmap(cbrfr,method=NULL, taxa.label = tl,taxa.order=tl)+ ylab("Relative abundance") + xlab(i) +scale_fill_gradient(name ="Relative abundance")
    plot(hm)
    # Saving the pdf
    dev.off()
    
    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2337, height = 1653, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    plot(hm)
    # Saving the pdf
    dev.off()
  }
  
  
}




############################    Phylum Level

#setwd("D:/Research/R codes")
#Variable aggregation
# different taxa ranks to iterate
taxal = c("Phylum")
# Looping through the selected sample variables
for (i in selectvs){
  # merge based on variables
  mergedphy =merge_samples(phyfiltered,group = i)
  
  # transform to relative abundance
  cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
  # Aggregation based on given taxa level
  for(tl in taxal){
    phylumGlommed = tax_glom(cbreedphy, tl, NArm=FALSE)
    
    # merging the OTU table with TAXA names and saving as a CSV
    # converting each table to dataframes
    ttdf=data.frame(as(tax_table(phylumGlommed), "matrix"))
    # have to transpose the OTU table as merging automatically transposes the matrix
    otdf=as.data.frame(t(otu_table(phylumGlommed)))
    # merging
    mtable= merge.data.frame(ttdf,otdf, by=0,all.x = F)
    # saving as a CSV
    csname= paste( tl, i,"relative_abandance.csv", sep="_")
    write.csv(mtable,csname)
    
    # filtering taxa based on percentage
    cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.01, TRUE)
    
    #Plotting the stacked barplots
    #File name
    trfname= paste( tl, i,"barplot.pdf", sep="_")
    jpegname= paste( tl, i,"barplot.jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    # Plotting the stacked barplots
    p = plot_bar(subset_taxa(cbrfr,!is.na(Phylum)),fill = tl) + ylab("Relative abundance") + xlab(i)
    p = p + geom_bar(color="black",stat="identity", position="stack")
    p=p+scale_fill_manual(values=c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
                                   "#9467bd","#8c564b","#e377c2","#7f7f7f",
                                   "#bcbd22","#17becf","#aec7e8","#ffbb78",
                                   "#98df8a","#ff9896","#c5b0d5","#c49c94"))
    plot(p)
    # Saving the pdf
    dev.off()
    
    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2337, height = 1653, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    plot(p)
    # Saving the pdf
    dev.off()
    
    
    
    # Heatmaps
    #File name
    trfname= paste( tl, i,"heatmap.pdf", sep="_")
    jpegname= paste( tl, i,"heatmap.jpeg", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    hm=plot_heatmap(cbrfr,method=NULL, taxa.label = tl,taxa.order=tl)+ ylab("Relative abundance") + xlab(i) + scale_fill_gradient(name ="Relative abundance")
    plot(hm)
    # Saving the pdf
    dev.off()

    # generate a high quality 350 ppi image
    jpeg(filename = jpegname,
         width = 2337, height = 1653, units = "px", pointsize = 12,
         quality = 75, bg = "white", res = 350)
    #plotting
    plot(hm)
    # Saving the pdf
    dev.off()
  }
  
  
}

#######################################           Pathogen Analysis          ######################
#############   pathogen - farm code bar plot
theme_set(theme_bw())
cbreedphy <- transform_sample_counts(phyfiltered, function(x) x / sum(x))

tax_table_data <- as.data.frame(cbreedphy@tax_table)
write.csv(tax_table_data, file = "tax_table.csv", row.names = FALSE)

gp.ch = subset_taxa(cbreedphy, Species == c("aerosaccus","agalactiae","nasimurium","saprophyticus"))

# saving as a CSV
csname= paste( "pathogens_farmcode.csv", sep="_")
write.csv(gp.ch,csname)

#File name
trfname= paste("RA_pathogen_FarmCode_barplot2.pdf")
jpegname= paste("RA_pathogen_FarmCode_barplot2.jpeg")
# Command to save as a pdf
pdf(file=trfname, paper="a4r",width = 0, height = 0)
# Plotting the stacked barplots
p=plot_bar(gp.ch,x="FarmCode", fill="Species")+ylab("Relative abundance")
p = p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot(p)
# Saving the pdf
dev.off()


# generate a high quality 350 ppi image
jpeg(filename = jpegname,
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot(p)
# Saving the pdf
dev.off()



##################   pathogen - sample bar plot
#setwd("D:/Research/R codes/impData")

theme_set(theme_bw())
cbreedphy <- transform_sample_counts(phyfiltered, function(x) x / sum(x))
gp.ch = subset_taxa(cbreedphy, Species == c("aerosaccus","agalactiae","nasimurium","saprophyticus"))

#File name
trfname= paste("RA_pathogen_Sample_barplot2.pdf")
jpegname= paste("RA_pathogen_Sample_barplot2.jpeg")
# Command to save as a pdf
pdf(file=trfname, paper="a4r",width = 0, height = 0)
# Plotting the stacked barplots
p=plot_bar(gp.ch,x="Sample", fill="Species")+ylab("Relative abundance")
p = p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot(p)

# Saving the pdf
dev.off()


# generate a high quality 350 ppi image
jpeg(filename = jpegname,
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot(p)
# Saving the pdf
dev.off()

########################################### Obtaining RA barplot cleanliness
#setwd("D:/Research/R codes/impData")

theme_set(theme_bw())
cbreedphy <- transform_sample_counts(phyfiltered, function(x) x / sum(x))
#gp.ch = subset_taxa(cbreedphy, Species == c("aerosaccus","agalactiae","nasimurium","saprophyticus"))
gp.ch = subset_taxa(cbreedphy, Species == c("aerosaccus","agalactiae","nasimurium","saprophyticus"))

#gp.ch = subset_taxa(cbreedphy, Species == c("saprophyticus"))


#File name
trfname= paste("RA_pathogen_cleanliness_barplot2.pdf")
jpegname= paste("RA_pathogen_cleanliness_barplot2.jpeg")
# Command to save as a pdf
pdf(file=trfname, paper="a4r",width = 0, height = 0)
# Plotting the stacked barplots
#p=plot_bar(gp.ch,x="FarmCode", fill="Species",facet_grid=~cleanliness)+ylab("Relative abundance")
p=plot_bar(gp.ch, fill="Species",facet_grid=~cleanliness)+ylab("Relative abundance")
p = p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot(p)

# Saving the pdf
dev.off()


# generate a high quality 350 ppi image
jpeg(filename = jpegname,
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot(p)
# Saving the pdf
dev.off()

#############################################################################
###########################                  Obtaining RA barplot FarmSize
#setwd("D:/Research/R codes/impData")

theme_set(theme_bw())
cbreedphy <- transform_sample_counts(phyfiltered, function(x) x / sum(x))
gp.ch = subset_taxa(cbreedphy, Species == c("aerosaccus","agalactiae","nasimurium","saprophyticus"))

#File name
trfname= paste("RA_pathogen_FarmSize_barplot2.pdf")
jpegname= paste("RA_pathogen_FarmSize_barplot2.jpeg")
# Command to save as a pdf
pdf(file=trfname, paper="a4r",width = 0, height = 0)
#p=plot_bar(gp.ch,x="FarmCode", fill="Species",facet_grid=~FarmSize)+ylab("Relative abundance")
p=plot_bar(gp.ch, fill="Species",facet_grid=~FarmSize)+ylab("Relative abundance")
p = p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot(p)

# Saving the pdf
dev.off()


# generate a high quality 350 ppi image
jpeg(filename = jpegname,
     width = 2337, height = 1653, units = "px", pointsize = 12,
     quality = 75, bg = "white", res = 350)
#plotting
plot(p)
# Saving the pdf
dev.off()
