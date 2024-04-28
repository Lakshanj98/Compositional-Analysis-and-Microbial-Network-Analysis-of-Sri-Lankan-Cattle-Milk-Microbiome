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
#Reading the meta data
metx=read.delim("metadata_table.txt",sep=",",header=T,row.names=sample_names(physeq))
head(metx)
#Meta data as sample data
metx=metx[,-1]
head(metx)
II
metx=sample_data(metx)
metx
#Merging the sample
complete_physeq=merge_phyloseq(physeq,metx)
complete_physeq
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
III
complete_physeq=ps_drop_incomplete(complete_physeq)
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
# filtering the taxa based on prevalence
# remove the phyla with total abundance(OTU count) less than 50
phyla2Filter = c("Synergistetes", "OP8", "TM6",
"SR1","GN02","Fibrobacteres",
"Nitrospirae","WPS-2","Planctomycetes")
IV
# Filter entries with unidentified Phylum.
# unknown <- subset_taxa(physeq = tax_table(complete_physeq), !is.na(Phylum) & !Phylum %in% c("", "NA"))
phyfiltered = subset_taxa(complete_physeq, !Phylum %in% phyla2Filter)
phyfiltered

############ Pathogens analysis
###########################################################################
# obtaining RA barplot cleanliness
#setwd("D:/Research/R codes/impData")
theme_set(theme_bw())
cbreedphy <- transform_sample_counts(phyfiltered, function(x) x / sum(x))
gp.ch = subset_taxa(cbreedphy, Species == c("aerosaccus","agalactiae","nasimurium","saprophyticus"))
#File name
trfname= paste("RA_pathogen_cleanliness_barplot2.pdf")
jpegname= paste("RA_pathogen_cleanliness_barplot2.jpeg")
# Command to save as a pdf
pdf(file=trfname, paper="a4r",width = 0, height = 0)
# Plotting the stacked barplots
p=plot_bar(gp.ch, x="FarmCode", fill="Species",facet_grid=~cleanliness)+ylab("Relative abundance")
p = p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot(p)
# Saving the pdf
dev.off()
# generate a high quality 350 ppi image
jpeg(filename = jpegname,
width = 2337, height = 1653, units = "px", pointsize = 12,
quality = 75, bg = "white", res = 350)
XXV
#plotting
plot(p)
# Saving the pdf
dev.off()
#############################################################################
# obtaining RA barplot FarmSize
#setwd("D:/Research/R codes/impData")
theme_set(theme_bw())
cbreedphy <- transform_sample_counts(phyfiltered, function(x) x / sum(x))
gp.ch = subset_taxa(cbreedphy, Species == c("aerosaccus","agalactiae","nasimurium","saprophyticus"))
#File name
trfname= paste("RA_pathogen_FarmSize_barplot2.pdf")
jpegname= paste("RA_pathogen_FarmSize_barplot2.jpeg")
# Command to save as a pdf
pdf(file=trfname, paper="a4r",width = 0, height = 0)
p=plot_bar(gp.ch, x="FarmCode", fill="Species",facet_grid=~FarmSize)+ylab("Relative abundance")
p = p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot(p)
# Saving the pdf
dev.off()
# generate a high quality 350 ppi image
jpeg(filename = jpegname,
width = 2337, height = 1653, units = "px", pointsize = 12,
XXVI
quality = 75, bg = "white", res = 350)
#plotting
plot(p)
# Saving the pdf
dev.off()
