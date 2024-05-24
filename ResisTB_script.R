# Load packages
library(vegan)
library(DESeq2)
library(qiime2R)
library(phyloseq)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(LTBI-ales)
library(data.table)
library(fBasics)
library(forcats)
library(dplyr)
library(ggLTBI-i)

# Load data
physeq<-qza_to_phyloseq("table-merge-filter.qza","rooted-tree_quality.qza","taxonomy.qza", "ResisTB_mapping_DMM3.txt")

# Remove taxa with 0 abundance
otu.table = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)

# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
otu.relative.table = transformSampleCounts(otu.table, normalizeSample)

# Check
rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))

# Subset
LTBI+.table = subset_samples(otu.relative.table, LTBI_status %in% c("LTBI+"))
LTBI-.table = subset_samples(otu.relative.table, LTBI_status %in% c("LTBI-"))
age.above.table = subset_samples(otu.relative.table, age_ %in% c("35-60"))
age.below.table = subset_samples(otu.relative.table, age_ %in% c("18-25"))


# Save
save.image(file="ResisTB.RData")

# Load
load(file="ResisTB.RData")

## ALPHA DIVERSITY

alpha <- estimate_richness(physeq)
write.csv(alpha, "alpha.diversity.csv")

#Figures made on Graphpad

## BETA DIVERSITY

## Figure 2B
# Create a distance matrix 
vegdist = distance(otu.relative.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdLTBI-ale <- cmdLTBI-ale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdLTBI-ale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdLTBI-ale, y = sample_data(otu.relative.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ LTBI_status,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="LTBI_status",suffixes=c("",".centroid"))

# Plot
pdf("Figure2B.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= LTBI_status)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  LTBI-ale_color_manual(values=c("forestgreen", "red")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= LTBI_status), size=0) +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= LTBI_status, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=LTBI_status, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics for Beta diversity
# Done as described below for all beta diversity comparisons including for variables LTBI_status,age_, Field_site, Current.CD4.count, INH.prophylaxis., Cotrimoxazole.prophylaxis.,BMI, DMM_pruned99.
  
bray = distance(, "bray")
adonis2(bray ~ sample_data(otu.relative.table)$LTBI_status)

## Figure 3B
# Create a distance matrix 
vegdist = distance(LTBI+.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdLTBI-ale <- cmdLTBI-ale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdLTBI-ale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdLTBI-ale, y = sample_data(LTBI+.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ age_,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="age_",suffixes=c("",".centroid"))

# Plot
pdf("Figure3B.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= age_)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  LTBI-ale_color_manual(values=c("mediumpurple1", "deeppink1")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= age_), size=0) +
  #ggtitle("Stool") +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= age_, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=age_, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Figure 3D
# Create a distance matrix 
vegdist = distance(LTBI-.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdLTBI-ale <- cmdLTBI-ale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdLTBI-ale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdLTBI-ale, y = sample_data(LTBI+.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ age_,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="age_",suffixes=c("",".centroid"))

# Plot
pdf("Figure3D.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= age_)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  LTBI-ale_color_manual(values=c("mediumpurple1", "deeppink1")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= age_), size=0) +
  #ggtitle("Stool") +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= age_, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=age_, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

## Figure 5B
# Create a distance matrix 
vegdist = distance(otu.relative.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdLTBI-ale <- cmdLTBI-ale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdLTBI-ale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdLTBI-ale, y = sample_data(otu.relative.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ DMM_pruned99,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="DMM_pruned99",suffixes=c("",".centroid"))

# Plot
pdf("Figure5B.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= DMM_pruned99)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  LTBI-ale_color_manual(values=c("blue", "orange", "purple")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= DMM_pruned99), size=0) +
  #ggtitle("Stool") +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= DMM_pruned99, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=DMM_pruned99, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

## DMM
#Figure S3
#load Libraries
library(DirichletMultinomial)
library(lattice)
library(xtable)
library(parallel)

#write data table
write.table(otu_table(otu.table),file="OTU.Table.txt",sep="\t", col.names = NA, row.names = TRUE)

#load data table
dmm.table = read.table('OTU.Table.txt', header=T, sep="\t", as.is=TRUE)

#convert table to matrix
count <- t(as.matrix(dmm.table, row.names=1, col.names=1))

#fix matrix with right colnames and convert to numeric from string
colnames(count) <- count[1,]
count<-count[-1,]
class(count) <-"numeric"

#Remove rows with 0 values for all data
count <- count[rowSums(count != 0)>0,]

#display matrix
count[1:5, 1:3]

#fit model for 1-7 clusters
fit.all <- mclapply(1:5, dmn, count=count, verbose=TRUE)

#save models to file
save(fit.all, file="fit.rda")
save.image(file="fit.RData")

lplc <- sapply(fit.all, laplace)

#plot the Lapplce based on number of clusters
pdf("FigureS3.pdf", height = 10, width = 10)
plot(lplc, type="b", xlab="Number of Dirichlet Components",
     ylab="Model Fit")
dev.off()

#save dataframe 
df <- as.data.frame(lplc)
df$fit <-1:nrow(df)

#display the number of clusters with the lowest laplace
(best <- fit.all[[lplc=3]])

#create table with grouping of clusters
grouping <- mixture(best, assign=TRUE)
write.table(grouping, file="grouping_DMM3.txt", sep="\t")

#save models to file
save(fit.all, file="fit.rda")
save.image(file="fit.RData")

##HEATMAP
#Figure5C
#set data tables
GenusData <-otu_table(Genus.Rel.table1B) #pruned to selected Genera based on abundance
#create vector to lable by Final_Sub_Analysis_Code_ca
colnames(sample_data(Genus.Rel.table1B))
sample_data(Genus.Rel.table1B)$DMM_pruned99
SampleVector = sample_data(Genus.Rel.table1B)$DMM_pruned99 #There are 2 different ones
#duplicate to create a color vector and replace value w/ color
#Colorvector can only replace numbers! ### Not true, can also take strings
Colorvector <-SampleVector
Colorvector <- replace(Colorvector, which (Colorvector == "C1"), "blue")
Colorvector <- replace(Colorvector, which (Colorvector == "C2"), "orange")
Colorvector <- replace(Colorvector, which (Colorvector == "C3"), "purple")


##Cluster Bray Heatmap
#cluster Genuses(row)
GenusData.Bray.dist <-vegdist(GenusData, method = "bray")
Genus.Bray.clus <-hclust(GenusData.Bray.dist, "aver")
Samples.Bray.dist = phyloseq::distance(GenusData, method="bray")
Samples.cluster.Bray = hclust(Samples.Bray.dist, "aver")
#Here we are able to change the names for genuses that are labelled as "g__" --> Come back to this
tax_table(Genus.Rel.table1B)
Genus.rel.table1B.New.Names = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.table1B))), ntaxa(Genus.Rel.table1B)), Genus.Rel.table1B)


# Add a new rank, Strain, with the Genus ids
tax_table(Genus.rel.table1B.New.Names) <- cbind(tax_table(Genus.rel.table1B.New.Names), Strain=taxa_names(Genus.rel.table1B.New.Names))
# Define the ranks you want to include
myranks = c("Class","Order", "Family", "Genus")
mylabels = apply(tax_table(Genus.rel.table1B.New.Names)[, myranks], 1, paste, sep="", collapse="_")
# Add concatenated labels as a new rank after strain
tax_table(Genus.rel.table1B.New.Names) <- cbind(tax_table(Genus.rel.table1B.New.Names), catglab=mylabels)
tax_table(Genus.rel.table1B.New.Names)
#Now Plot Heat map with dendograms
mypalette <- colorRampPalette(c('#FFFFFF','#4169E1','#0000CD'))
pdf("Figure5C.pdf", height = 10, width = 10)
heatmap.2(GenusData,
          density.info = "none",
          trace = "none",
          dendrogram = "both",
          Rowv = as.dendrogram(Genus.Bray.clus),
          Colv = as.dendrogram(Samples.cluster.Bray),
          labRow=tax_table(Genus.rel.table1B.New.Names)[,"catglab"],
          cexRow = .6,
          labCol = sample_data(Genus.Rel.table1B)$SampleID,
          cexCol = .8,
          col = mypalette(17),
          symm=F,symkey=F,symbreaks=T, LTBI-ale="none",
          breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
          ColSideColors= Colorvector,
          main = "Heatmap of Bray-Curtis Distance",
)
dev.off()

##DIFFERENTIAL ABUNDANCE

library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(LTBI-ales)
library(data.table)
library(fBasics)
library(forcats)
library(vegan)
library(dplyr)
library(MetaboSignal)
library(phyloseq)

theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
              legend.position="none")

# Set cut-off
alpha = 0.20
alpha <- 0.20

#Figure 2C
#subset otu table to genus level
subset.genus.table = tax_glom(otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)


# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ LTBI_status)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$LTBI_status<- droplevels(diagdds$LTBI_status)

# choose variable of interest and reference
diagdds$LTBI_status<- relevel(diagdds$LTBI_status, ref ="LTBI-")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, LTBI_status%in% c("LTBI+"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, LTBI_status%in% c("LTBI-"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_LTBI_status.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "forestgreen"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure2C.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure 3C

#subset otu table to genus level
subset.genus.table = tax_glom(LTBI+.otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI+.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI+.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)


# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ age_)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$age_<- droplevels(diagdds$age_)

# choose variable of interest and reference
diagdds$age_<- relevel(diagdds$age_, ref ="18-25")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, age_%in% c("35-60"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, age_%in% c("18-25"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_LTBI+_age.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "deeppink1"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumpurple1"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure3C.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure 3F

#subset otu table to genus level
subset.genus.table = tax_glom(LTBI-.otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI-.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI-.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)


# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ age_)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$age_<- droplevels(diagdds$age_)

# choose variable of interest and reference
diagdds$age_<- relevel(diagdds$age_, ref ="18-25")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, age_%in% c("35-60"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, age_%in% c("18-25"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_LTBI-_age.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "deeppink1"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumpurple1"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure3F.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure 5D

subset.otu.table = subset_samples(otu.table, DMM_pruned99 %in% c("C1", "C2"))
nsamples(subset.otu.table)

#subset otu table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)


# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ DMM_pruned99)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$DMM_pruned99<- droplevels(diagdds$DMM_pruned99)

# choose variable of interest and reference
diagdds$DMM_pruned99<- relevel(diagdds$DMM_pruned99, ref ="C1")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, DMM_pruned99%in% c("C2"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, DMM_pruned99%in% c("C1"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_C1vLTBI-2.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "orange"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "blue"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure5D.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure 5E

subset.otu.table = subset_samples(otu.table, DMM_pruned99 %in% c("C1", "C2"))
nsamples(subset.otu.table)

#subset otu table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ DMM_pruned99)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$DMM_pruned99<- droplevels(diagdds$DMM_pruned99)

# choose variable of interest and reference
diagdds$DMM_pruned99<- relevel(diagdds$DMM_pruned99, ref ="C1")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, DMM_pruned99%in% c("C2"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, DMM_pruned99%in% c("C1"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_C1vLTBI-2.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "purple"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "blue"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure5E.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure 5F

subset.otu.table = subset_samples(otu.table, DMM_pruned99 %in% c("C2", "C3"))
nsamples(subset.otu.table)

#subset otu table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ DMM_pruned99)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$DMM_pruned99<- droplevels(diagdds$DMM_pruned99)

# choose variable of interest and reference
diagdds$DMM_pruned99<- relevel(diagdds$DMM_pruned99, ref ="C3")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, DMM_pruned99%in% c("C2"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, DMM_pruned99%in% c("C3"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_C3vLTBI-2.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "orange"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "purple"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure5F.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure 6A

#subset otu table to genus level
subset.genus.table = tax_glom(LTBI+.otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ IGRA_group)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$IGRA_group <- droplevels(diagdds$IGRA_group)

# choose variable of interest and reference
diagdds$IGRA_group<- relevel(diagdds$IGRA_group, ref ="below")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, IGRA_group%in% c("above"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, IGRA_group%in% c("below"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_LTBI+.IGRA_group.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "goldenrod"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "grey"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure6A.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure 6B

#subset otu table to genus level
subset.genus.table = tax_glom(LTBI+.otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ TST_group)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$TST_group <- droplevels(diagdds$TST_group)

# choose variable of interest and reference
diagdds$TST_group<- relevel(diagdds$TST_group, ref ="below")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, TST_group%in% c("above"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, TST_group%in% c("below"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_TST_group LTBI+.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "goldenrod"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "grey"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure6B.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure S1

#subset otu table to genus level
subset.genus.table = tax_glom(otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ INH.prophylaxis.)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$INH.prophylaxis. <- droplevels(diagdds$INH.prophylaxis.)

# choose variable of interest and reference
diagdds$INH.prophylaxis.<- relevel(diagdds$INH.prophylaxis., ref ="No")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, INH.prophylaxis.%in% c("Yes"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, INH.prophylaxis.%in% c("No"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_INH.prophylaxis.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "saddlebrown"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "salmon"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="FigureS1.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure S5A

#subset otu table to genus level
subset.genus.table = tax_glom(otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ IGRA_group)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$IGRA_group <- droplevels(diagdds$IGRA_group)

# choose variable of interest and reference
diagdds$IGRA_group<- relevel(diagdds$IGRA_group, ref ="below")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, IGRA_group%in% c("above"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, IGRA_group%in% c("below"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_IGRA_group_overall.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "goldenrod"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "grey"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="FigureS5A.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

#Figure S5B

#subset otu table to genus level
subset.genus.table = tax_glom(otu.table, taxrank = "Genus")
ntaxa(subset.genus.table)
#Prune data to less than 100 genera remaining
LTBI_status.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.10 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(LTBI_status.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

#Normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Change variable of interest
# you can adjust analyis by other variables e.g. Put confounders first
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ TST_group)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
1
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# change variable of interest
diagdds$TST_group <- droplevels(diagdds$TST_group)

# choose variable of interest and reference
diagdds$TST_group<- relevel(diagdds$TST_group, ref ="below")

diagdds<- DESeq(diagdds)
res <- results(diagdds)
res = res[order(res$padj, na.last = NA), ]
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,rownames(res))
res$row2 <- gsub('\\s+', '|', res$row2)
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Kingdom","Phylum","Class","Order","Family","Genus","Species","row2"))
res$names <- res$Taxa
res$Taxa <- res$row2
otu.to.save <-as.character(res$names)

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, TST_group%in% c("above"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, TST_group%in% c("below"))

experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)

control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

res$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Change abundance column headers
res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# change table name
write.table(res.1,file="abundance_0.20_TST_group_overall.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res$sig <- -log10(res$adj.P.Val)

sum(is.infinite(res$sig))

cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"

# Change colours
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "goldenrod"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "grey"

res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="FigureS5B.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 50 * res$abundance.experiment, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 50 * res$abundance.control,2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.20,max.overlaps = 50) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

##PICRUST
#Figure 4
#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(LTBI-ales)
library(data.table)
library(fBasics)
library(forcats)
library(vegan)
library(dplyr)
library(MetaboSignal)
library(phyloseq)
require(ade4)
require(vegan)
require(tidyverse)
require(adegenet, quietly = T)
require(plotly)
require(lubridate)
```

kegg_brite_map <- read.table("picrust1_KO_BRITE_map.tsv",
                              header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)

 
### When reading in tab-delimited file of KO predictions (PICRUSt2 output):
test_ko <- read.table("path_abun_unstrat_deLTBI-rip.txt", header=TRUE, sep="\t", row.names=1)
write.table(test_ko_L3,file="test_ko_L3.csv",sep=",", col.names = NA, row.names = TRUE)

#load counts
count_data = read.csv(file = "test_ko_L3.csv", header = T, sep = ",", row.names=1)
head(count_data)

# Convert to matrix
cts <- as.matrix(cts)
head(cts)

any(is.na(cts))
all(is.numeric(cts))

#Load metadata
sampleInfo <- read.table("ResisTB_mapping_DMM3.txt", header=T, sep="\t", row.names=1)
sampleInfo <- as.matrix((sampleInfo))
colnames(sampleInfo)
rownames(sampleInfo)

#set your alpha
alpha <- 0.20

#generate the DESeqDataSet - the core object of DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sampleInfo,
                              design = ~ LTBI_status)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(assay(dds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = estimateDispersions(dds)

#Drop rows with no data in your comparison variable
dds$LTBI_status <- droplevels(dds$LTBI_status)

#Choose which is the 'control' in your comparison variable --> positive is upregulated in LTBI+, negative is down-regulated
dds$LTBI_status <- relevel(dds$LTBI_status, ref ="LTBI-")

#Run the differential Analysis
dds<- DESeq(dds)

#Output the results from DESEQ into a table
res <- results(dds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

######get abundance data 
#decide what IDs to save 

#otu.to.save <-as.character(rownames(res4))
otu.to.save <-as.character(rownames(res))

#from main table we should get the mean across the row of the table
ko.table.df <- data.frame(assay(dds))
ko.table.df.meanRA <- rowMeans(ko.table.df)

#need to subset AND reorder just the IDs that we have
ko.table.df.meanRA.save <- ko.table.df.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- ko.table.df.meanRA.save

#subset relative table for groups being analysed to get their separate abundances
## Extract and subset count data
countdata = assay(dds)
coldata = colData(dds)
LTBI+.pruned.genus.rel.table = cts[, coldata$LTBI_status %in% c("LTBI+")]
LTBI-.pruned.genus.rel.table = countdata[, coldata$LTBI_status %in% c("LTBI-")]

#from relative table we should get the mean across the row of the otu table
LTBI+.pruned.genus.rel.table.df <- data.frame(LTBI+.pruned.genus.rel.table)
LTBI+.pruned.genus.rel.table.df.meanRA <- rowMeans(LTBI+.pruned.genus.rel.table.df)

LTBI-.pruned.genus.rel.table.df <- data.frame(LTBI-.pruned.genus.rel.table)
LTBI-.pruned.genus.rel.table.df.meanRA <- rowMeans(LTBI-.pruned.genus.rel.table.df)

#need to subset AND reorder just the otus that we have 
LTBI+.pruned.genus.rel.table.df.meanRA.save <- LTBI+.pruned.genus.rel.table.df.meanRA[otu.to.save]
LTBI-.pruned.genus.rel.table.df.meanRA.save <- LTBI-.pruned.genus.rel.table.df.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.LTBI+ <- LTBI+.pruned.genus.rel.table.df.meanRA.save
res$abundance.LTBI- <- LTBI-.pruned.genus.rel.table.df.meanRA.save

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Pathway.key","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "abundance", "abundance.dTBLs", "abundance.nTBLs"))
write.table(res,file="res__PICRUSt2.txt", sep="\t", col.names = NA, row.names = TRUE)

# Compute FDR in a log LTBI-ales
res$sig <- -log10(res$adj.P.Val)

#See how many are now infinite
sum(is.infinite(res$sig))

####If there infinite value set a maximum value for adj p value for the graph (e.g. 350)
res[is.infinite(res$sig),"sig"] <- 350

##Set the colors for your volcano plat
cols <- denLTBI-ols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "forestgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

pdf(file="Figure4.pdf", width=10, height=10)
ggplot(res, aes(x = logFC, y = sig,label=Pathway.key)) +
  geom_point(color=cols, size = ifelse(res$logFC>=0 & res$adj.P.Val < alpha,  2000 * res$abundance.LTBI+, ifelse(res$logFC<=-0 & res$adj.P.Val < alpha, 2000 * res$abundance.LTBI-,2)),alpha=0.7) + #Chose Colors and size for dots
#  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha, as.character(res$Pathway.key),'')),size=3,force=25, segment.colour="black",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC>0.5, as.character(res$Pathway.key),'')),size=2,force=25,segment.colour="black",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC<=-0.5, as.character(res$Pathway.key),'')),size=2,force=25,segment.colour="black",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

# Figure S2
# Bdivs distances were plotted on graphpad.

## LEFSE
#Figure S4
#load packages
library(microbiomeMarker)
library(purrr)

genus.table = tax_glom(otu.table, taxrank = "Genus")

#normalize
normalizeSample = function(x) {
  x/sum(x)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# Load/Save/wd-----


# set threshold ----
pvalue.thres <- 0.2
lda.thres <- 2

ps <- genus.table
colnames(tax_table(ps))=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(ps))

ps.rel <- transformSampleCounts(ps, normalizeSample)

#### Using wrt-----
set.seed(123) 
res <- run_lefse(ps,
                 group="DMM_pruned99",
                 subgroup = "none",
                 taxa_rank = "none", # this will give ASV values. Cannot use this with Cladogram
                 # needs to be "all" to run Cladogram (needs names listed as f__, g__, s__)
                 # needs to be "all"none" to do Bubble Plot (needs ASV) -- will somehow need to get ASV changed to names
                 transform = "identity", #log10, log10p
                 norm = "CPM", #CPM is method of choice for normalizing data for LEFSE analysis
                 kw_cutoff = 0.2,
                 #wilcoxon_cutoff = 0.05,
                 lda_cutoff = 2,
                 bootstrap_n = 30, #30 is default
                 bootstrap_fraction = 2/3,
                 multigrp_strat = TRUE
)
res1 <- data.frame(marker_table(res))
res1$marker <- rownames(res1)
res1$number <- c(1:nrow(res1))

taxtable <-data.frame(tax_table(ps))

taxtable$names <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",taxtable$Species,sep=""),
                               ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,sep=""),
                                      ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",sep=""),
                                             ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",sep=""),
                                                    ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",sep=""),
                                                           ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",sep=""),
                                                                  ifelse(!is.na(taxtable$Domain), paste("k_",taxtable$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(taxtable))))))))))


#creates column "nameASV" that includes trailing ASV
taxtable$nameASV <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",rownames(taxtable),sep=""),
                                 ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,"_",rownames(taxtable),sep=""),
                                        ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",rownames(taxtable),sep=""),
                                               ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",rownames(taxtable),sep=""),
                                                      ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",rownames(taxtable),sep=""),
                                                             ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",rownames(taxtable),sep=""),
                                                                    ifelse(!is.na(taxtable$Domain), paste("k_",taxtable$Domain, "_p__c__o__f__g__",rownames(taxtable),sep=""), paste(rownames(taxtable))))))))))

#create column ASV that is from row name
taxtable$asv <- rownames(taxtable)


merged <- merge(x=res1, y=taxtable, by.x="feature",by.y="asv", all.x=T)

merged <- merged[order(merged$number),]
res1$asv <- res1$feature
res1$feature <- merged$names
res1$number <- merged$number


res1 <- res1[order(res1$enrich_group, res1$ef_lda),]
res1$feature <- factor(res1$feature, levels = res1$feature)

#Subset LEFSE Table based on each comparison
resc1 <- res1[res1$enrich_group=="C1",]
resc2 <- res1[res1$enrich_group=="C2",] 
resc3 <- res1[res1$enrich_group=="C3",] 

#Subset OTU Table based on each comparison
otuC1 <- subset_samples(ps, DMM_pruned99 %in% c("C1"))
otuC2 <- subset_samples(ps, DMM_pruned99 %in% c("C2"))
otuC3 <- subset_samples(ps, DMM_pruned99 %in% c("C3"))

#Relative Abundance of subjset OTU Tables
## To normalize data you need to set a function
rel.otuC1 = transformSampleCounts(otuC1, normalizeSample)
rel.otuC2 = transformSampleCounts(otuC2, normalizeSample)
rel.otuC3 = transformSampleCounts(otuC3, normalizeSample)
#decide what otu to save based on LEFSE
otu.to.save1 <-as.character(resc1$asv)
otu.to.save2 <-as.character(resc2$asv)
otu.to.save3 <-as.character(resc3$asv)

#convert Relative Abundance Table to a dataframe
df.1.df <- data.frame(otu_table(rel.otuC1))
df.2.df <- data.frame(otu_table(rel.otuC2))
df.3.df <- data.frame(otu_table(rel.otuC3))
#from relative Abundance tables we should get the mean across the row of the otu table
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)
df.3.meanRA <- rowMeans(df.3.df)

#need to subset AND reorder just the otus that we are interested in from LEFSE 
df.1.meanRA.save <- df.1.meanRA[otu.to.save1]
df.2.meanRA.save <- df.2.meanRA[otu.to.save2]
df.3.meanRA.save <- df.2.meanRA[otu.to.save3]

#add the abundnace data for the res dataframe
resc1$abundance <- df.1.meanRA.save
resc2$abundance <- df.2.meanRA.save
resc3$abundance <- df.3.meanRA.save

#Merge the comparisons into one table
RESDF <- rbind(resc1,resc2,resc3)
#Calculate FDR & Subset for set FDR
RESDF$padjust <- p.adjust(RESDF$pvalue,method="BH")
RESDF <- data.frame(merge(x = RESDF, y = res1))

#Remove any NAs
RESDF <- RESDF[!is.na(RESDF$abundance),]


RESDF <-RESDF[RESDF$padjust < pvalue.thres & RESDF$ef_lda >= lda.thres, , drop=TRUE]
# need to order RESDF here so that RESDF$contamcolor works ...
RESDF <- RESDF[with
               (RESDF, order(enrich_group, ef_lda)),]


ggplot_LEFSE_bar <- function(obj) {
  # if (any(obj$padjust < pvalue.thres & any(obj$ef_lda >= lda.thres)))
  # obj[obj$padjust < pvalue.thres & obj$ef_lda >= lda.thres, , drop=TRUE] %>%
  obj %>%
    ggplot(., aes(x = ef_lda, y = feature, fill=enrich_group), alpha = 0.8)+
    geom_bar(stat="identity")+
    scale_fill_manual(name = "Clusters",
                      values=c("C1"="blue","C2"="goldenrod", "C3"="purple"),
                      labels=c("C1","C2", "C3"))+
    guides(fill= guide_legend(reverse=TRUE))+
    theme_pubr()+
    theme(
      axis.text = element_text(face="bold"),
      plot.title = element_text(face="bold"),
      axis.title = element_text(face="bold"),
      axis.text.y = element_text(angle=0,hjust=1, face = "bold.italic"))+
    xlab("LDA (log10)")+
    xlim(0,ceiling(max(obj$ef_lda)))+
    ylab("")
  # +
  # ggtitle(paste("LEfSe \nK-W FDR<",pvalue.thres,", LDA>=",lda.thres)) 
}

LEFSE_bar_genus <- ggplot_LEFSE_bar(RESDF)
LEFSE_bar_genus


# creating abundance plots
RESDF <- as_tibble(RESDF)
RESDF$padj <- as.numeric(RESDF$padj)
RESDF$padjust <- as.numeric(RESDF$padjust)
RESDF$pvalue <- as.numeric(RESDF$pvalue)
RESDF$ef_lda <- as.numeric(RESDF$ef_lda)
RESDF$abundance <- as.numeric(RESDF$abundance)
RESDF$enrich_group <- factor(RESDF$enrich_group, levels=c("C1","C2","C3"))
RESDF$feature <- factor(RESDF$feature, levels= RESDF$feature)

pvalue.thres
lda.thres

RESDF <- RESDF[RESDF$padjust < pvalue.thres & RESDF$ef_lda >= lda.thres,] # keep rows with pvalues < pvalue.thres

#order RESDF
RESDF <- RESDF[with
               (RESDF, order(enrich_group, ef_lda)),]


#order names in RESDF
RESDF$feature <- factor(RESDF$feature, levels=RESDF$feature)


abund <- data.frame(otu_table(ps.rel))
abund.select <- subset(abund, row.names(abund) %in% RESDF$asv)
RESDF <- subset(RESDF, RESDF$asv %in% row.names(abund.select))

dmm_pruned99 <- data.frame(sample_data(ps.rel))
dmm_pruned99 <- subset(dmm_pruned99, select = c("DMM_pruned99"))
dmm_pruned99$DMM_pruned99 <- factor(dmm_pruned99$DMM_pruned99, levels = c("C1","C2","C3"))

abund.select.transpose <- t(abund.select)
rownames(abund.select.transpose) <- gsub("X", "", rownames(abund.select.transpose))

# Replace . with - in all column names of abund.select
#for (i in 1:length(row.names(abund.select.transpose))) {
  #row.names(abund.select.transpose)[i] <- gsub(".", "-", row.names(abund.select.transpose)[i], fixed = TRUE)
}

# View the updated column names of abund.select
row.names(abund.select.transpose)


test <- merge(x = abund.select.transpose, y = dmm_pruned99, by = 0, all.x = T)

test1 <- data.frame(t(test))
row.names <- test1[,1]
names(test1) <- test1[nrow(test1),]
test1 <- test1[-1,]
test1 <- test1[-nrow(test1),]
test1 <- data.frame(sapply(test1[,1:ncol(test1)],as.numeric))
row.names(test1) <- colnames(test[,2:(ncol(test)-1)])
test1$feature <- row.names(test1)
test1 <- merge(x = test1, y = RESDF, by.x="feature", by.y = "asv", all.x = TRUE)
row.names(test1) <- test1$names

test1 <- test1[
  with(test1, order(enrich_group, ef_lda)),]

names(test1)
test1 <- subset(test1, select=-c(feature, enrich_group, ef_lda, pvalue, padj, marker, number,abundance, padjust))

test1 <- data.frame(t(test1))
names(test1) <- test1[nrow(test1),]
test1 <- test1[-nrow(test1),]

test1$DMM_pruned99<- row.names(test1)
test1$DMM_pruned99[test1$DMM_pruned99 %in% paste0("C1", )] <- "C2"<- "C3"
test1$DMM_pruned99 <- gsub("C1", "C2", "C3" test1$DMM_pruned99)

test1[grepl("C1",test1$DMM_pruned99),]$DMM_prune
test1[grepl("Cases Uninvolved BALF",test1$DMM_pruned99),]$DMM_pruned99<- "Cases Uninvolved BALF"
test1$DMM_pruned99 <- factor(test1$DMM_pruned99, levels=c("C1","C2","C3"))
test1[,1:(ncol(test1)-1)] <- data.frame(sapply(test1[,1:(ncol(test1)-1)],as.numeric))

##### mean abundance plot ----
a_genus <- tidyr::gather(test1, key="id",value="value",1:(ncol(test1)-1), factor_key=TRUE)
a_genus <- merge(x = a_genus, y = RESDF, by.x = "id", by.y = "feature", all.x = T)
a_genus <- subset(a_genus, select = -c(enrich_group, ef_lda, pvalue, padj, marker, number, asv, abundance, padjust))
a_genus <- as_tibble(a_genus)
a_genus


a_genus <- a_genus %>%
  dplyr::group_by(id, DMM_pruned99) %>%
  mutate(count = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = sd(value)) %>%
  mutate(se = sd/sqrt(n())) %>%
  mutate(ic=se * qt((1-0.05)/2 + .5, count-1 )) %>%
  mutate(meanpercent = mean * 100) %>%
  ungroup() %>%
  distinct(.,id,DMM_pruned99,count,mean,sd,se,ic, .keep_all= TRUE)
a_genus

p_lefse_relab_mean_genus_log <- a_genus %>%
  ggplot(., aes(x = meanpercent, y = id, fill = DMM_pruned99)) + 
  geom_bar(stat="identity", position=position_dodge(0.9), width = 0.9)+
  geom_bar(stat="identity", position=position_dodge(0.9), width = 0.9, color="black", show_guide=FALSE, size = 0.2)+
  geom_errorbar( aes(xmin=(mean-se)*100, xmax=(mean+se)*100, group = DMM_pruned99), width=0.4, colour="black", alpha=0.9, size=.25, position = position_dodge(0.9))+
  scale_fill_manual(name = "Clusters",
                    values=c("C1"="blue","C2"="goldenrod","C3"="Purple"),
                    labels=c("C1","C2", "C3"))+
  guides(fill= guide_legend(reverse=TRUE))+
  theme_classic()+
  theme(axis.text.y = element_text(angle=0,hjust=1))+
  xlab("Relative Abundance (%)")+
  ylab("")+
  ggtitle("")+
  scale_x_continuous(trans= scales::pseudo_log_trans(base=10),
                     breaks = c(0, 5, 20, 50, 100),
                     labels = c(0, 5, 20, 50, 100),
                     limits = c(0, 100),
                     expand = c(0, 0)
  )
p_lefse_relab_mean_genus_log

pdf("LEFSE.pdf", height = 4, width=8)
LEFSE_bar_abundance_mean_log_genus<- ggarrange(LEFSE_bar_genus + geom_vline(xintercept=lda.thres,linetype="longdash",color="black",size =0.5),
                                               p_lefse_relab_mean_genus_log + 
                                                 theme(axis.text.y=element_blank(),
                                                       axis.title=element_text(face="bold"),
                                                       axis.text.x=element_text(colour="black", face="bold"),
                                                       axis.ticks=element_line(colour="black"),
                                                       legend.background = element_rect(color=NA),
                                                       legend.position = "top"),
                                               ncol=2,
                                               nrow=1,
                                               align = "h",
                                               labels =c("A","B"),
                                               common.legend=TRUE,
                                               legend = "bottom",
                                               widths = c(2,1)
)
LEFSE_bar_abundance_mean_log_genus
ggsave("LEFSE.pdf", LEFSE_bar_abundance_mean_log_genus)
dev.off()