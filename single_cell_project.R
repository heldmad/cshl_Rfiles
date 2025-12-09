# Chapter 1. Basic scRNA-seq workflow in Monocle3, some ggplot, and why you must 
# explore your data! 
# Mary B. O'Neill & Anh Vo, BBI BAT-Lab, Single Cell Genomics
# CSHL Computational Genomics Course 2025

################################################################################
###### Setup ###################################################################
################################################################################
# Set your working directory to where the data is downloaded
setwd("~/cshl_Rfiles/cshl_Rfiles/") #CHANGE ME!!!

library(monocle3)
library(ggplot2)
library(tidyverse)
library(viridis)
library(data.table)
library(scuttle)
library(cowplot)
library(edgeR)
library(scales)
library(DropletUtils)
library(SoupX)

################################################################################
###### Let's start with a small, publicly available dataset ####################
###### from 10X Genomics website; 1K PBMCs #####################################
################################################################################
# Read the data
#lung <- load_cellranger_data("outs")
#?load_cellranger_data
#?load_mm_data

#how to load in common gene names from their ensembl code: go to ensembl, go to 
#MOUSE GENES dataset, select features c("Gene Stable ID", "Gene Name", 
#etc. anything else you need)

## An alternative method of loading the data
lungset <- load_mm_data(mat_path = "raw_feature_bc_matrix/matrix.mtx.gz", 
                        feature_anno_path = "raw_feature_bc_matrix/features.tsv.gz", 
                        cell_anno_path = "raw_feature_bc_matrix/barcodes.tsv.gz")

# learn about the cell_data_set class
class(lungset) #cell_data_set
?cell_data_set
lungset #dim: 32285 7788 

# multiple ways to access the column data
head(pData(lungset))
head(colData(lungset))

# multiple ways to access the row data
head(fData(lungset))
head(rowData(lungset))

#change column names to have one column named gene_short_name
colnames(rowData(lungset)) <- c("gene_short_name", "type")
head(rowData(lungset))
# multiple ways to access the assays
counts(lungset)[1:10,1:10]
assay(lungset, 'counts')[1:10,1:10]

# Pre-process the data
?preprocess_cds
lungset <- preprocess_cds(lungset, method = 'PCA', num_dim = 50)
lungset #reducedDimNames(1): PCA

# Plot PC variance
plot_pc_variance_explained(lungset)

# Reduce dimensionality - default params
?reduce_dimension
lungset <- reduce_dimension(lungset)
lungset #reducedDimNames(2): PCA UMAP

# Reduce dimensionality - minDist = 0.05, nn = 5
#?reduce_dimension
#lungset_nn5 <- reduce_dimension(lungset, 
#                                umap.min_dist = 0.01,
#                                umap.n_neighbors = 5)
#lungset_nn5 #reducedDimNames(2): PCA UMAP

#save in object for easier plotting later
lungset$UMAP1 <- reducedDim(lungset, "UMAP")[,1]
lungset$UMAP2 <- reducedDim(lungset, "UMAP")[,2]

#lungset_nn5$UMAP1 <- reducedDim(lungset_nn5, "UMAP")[,1]
#lungset_nn5$UMAP2 <- reducedDim(lungset_nn5, "UMAP")[,2]

# Plot within Monocle3
?plot_cells
plot_cells(lungset)
plot_cells(lungset, genes= c("Wt1", "Acta2", "Malat1"))


# Run an unsupervised clustering
?cluster_cells
lungset <- cluster_cells(lungset)
lungset$clusters <- clusters(lungset) #save this in colData for easier access
plot_cells(lungset)


# Marker gene detection - data driven
?top_markers
marker_test_res <- top_markers(lungset, group_cells_by="cluster")

?saveRDS
saveRDS(marker_test_res, file="lungset_marker_test_results.RDS") #save for later use

str(marker_test_res) #get in the habit of looking at data structures
head(marker_test_res) #show the first 6 lines

# filter the list to retain 3 top markers for each group
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

# reference the gene_id
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

# make a dot plot of data-driven top markers and some cannonical marker genes
?plot_genes_by_group
plot_genes_by_group(lungset,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

################################################################################
###### Quality Control #########################################################
################################################################################

# Quality Control (QC) is a trade off. Being too aggressive can results in the 
# loss of rare cell populations, being too permissive can complicate cell 
# annotation and differential expression analyses. EVERY DATASET IS UNIQUE!

# The heart of scRNA-seq data is the feature-by-cell matrix
lungset
counts(lungset)[1:10,1:10] #sparse matrix, "." = 0; single-cell data is zero inflated! 

?detect_genes
lungset <- detect_genes(lungset)
head(fData(lungset))
head(pData(lungset))

expressed <- data.frame(rowData(lungset)) %>% arrange(desc(num_cells_expressed))
counts(lungset)[head(expressed$id, n=10), 1:10]

# Calculate mito content
head(fData(lungset))

# How might we identify mitochondrial genes?
fData(lungset)$MT <- grepl("^mt-", rowData(lungset)$gene_short_name)
table(fData(lungset)$MT) #sanity check

pData(lungset)$MT_reads <- Matrix::colSums(exprs(lungset)[fData(lungset)$MT,])
pData(lungset)$MT_perc <- pData(lungset)$MT_reads/Matrix::colSums(exprs(lungset))*100
summary(lungset$MT_perc)


pData(lungset)$n.umi <- Matrix::colSums(exprs(lungset))

head(pData(lungset))

# Lets examine some common QC metrics & learn about ggplot, a very flexible
# R package for plotting

# Your data must be a dataframe if you want to plot it using ggplot2
# We achieve this by forcing the column data to a dataframe
# The aes function specifies how variables in your dataframe map to features on your plot
# Geoms specify how your data should be represented 

# Most basic
ggplot(data.frame(colData(lungset)), aes(x=n.umi, y=num_genes_expressed, color=MT_perc)) +
  geom_point()

# Let's spruce this up a bit
ggplot(data.frame(colData(lungset)), aes(x=n.umi, y=num_genes_expressed, color=MT_perc)) +
  geom_point() +
  theme_light() + #there are several preset themes, you can also customize 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + #fancy log scale
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) #fancy log scale

# And color by cluster  
ggplot(data.frame(colData(lungset)), aes(x=n.umi, y=num_genes_expressed, color=clusters)) +
  geom_point() +
  theme_light() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))


# We can use other types of plots (and layer!): violin, boxplot
ggplot(data.frame(colData(lungset)), aes(x=clusters, y=n.umi, fill=clusters)) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_light()

summary(lungset$n.umi)

# And the % mitochondrial
ggplot(data.frame(colData(lungset)), aes(y=MT_perc)) +
  geom_density(fill="salmon") + #yet another type of plot
  coord_flip() +
  theme_light()

summary(lungset$MT_perc)

# Let's ensure there are no cells with high mito & high umi
ggplot(data.frame(colData(lungset)), aes(x=n.umi, y=MT_perc)) +
  geom_point(alpha=0.3) + 
  theme_light() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

# This is usually the case, but sometimes a certain cell type or cell state
# is associated with high mitochondrial content in a biologically meaningful
# context and we do not want to throw those cells out.


#random QC methods, looking at clustering of cells and coloring based on MT percentage and n.umis and n.genes expressed
hist(colData(lungset)$n.umi)
hist(colData(lungset)$num_genes_expressed)

plot_cells(lungset, color_cells_by = "n.umi")
plot_cells(lungset, color_cells_by = "num_genes_expressed")
plot_cells(lungset, color_cells_by = "MT_perc")

lungset$qc <- ifelse(lungset$n.umi > 500 & lungset$MT_perc < 20, "KEEP", "REMOVE")

plot_cells(lungset, color_cells_by = "qc")
plot_cells(lungset, color_cells_by = "clusters")

lungsetqc <- lungset[, lungset$qc == "KEEP"]

cluster29 <- lungset[]


# What thresholds should be used? One approach is to use adaptive thresholds.
# If you can reasonably assume that most cells are of acceptable quality,
# identifying and removing outliers may be a solid approach.

?isOutlier


# The number of UMIs are not normally distributed so we use the log. Play with
# the number of MADs and see how this changes.
lib <- isOutlier(lungset$n.umi, 
                 log=TRUE, 
                 nmads=3, #change me
                 type=c("both"))
attr(lib, "thresholds") #128 UMIs at 3 MADs, 90 at 5, 182 at 1

genes <- isOutlier(lungset$num_genes_expressed, 
                   log=TRUE,
                   nmads=3, #change me
                   type=c("both"))
attr(genes, "thresholds") #103 genes at 3 MADs, 72 at 5, 146 at 1

high.mito <- isOutlier(lungset$MT_perc, 
                       nmads=5, #change me
                       type=c("higher"))
attr(high.mito, "thresholds") #2% at 3, 3% at 5, 1% at 5

# Let's plot those thresholds!
#by mito 
ggplot(data.frame(colData(lungset)), aes(x=n.umi, y=num_genes_expressed, color=MT_perc)) + #also try coloring by clusters
  geom_point(alpha=0.3) +
  theme_light() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept=attr(lib, "thresholds")[1], linetype="dotted", color="red") +
  geom_vline(xintercept=attr(lib, "thresholds")[2], linetype="dotted", color="red") +
  geom_hline(yintercept=attr(genes, "thresholds")[1], linetype="dotted", color="blue") +
  geom_hline(yintercept=attr(genes, "thresholds")[2], linetype="dotted", color="blue")

#by cluster - should be 30 colors
ggplot(data.frame(colData(lungset)), aes(x=n.umi, y=num_genes_expressed, color=clusters)) + #also try coloring by clusters
  geom_point(alpha=0.3) +
  theme_light() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept=attr(lib, "thresholds")[1], linetype="dotted", color="red") +
  geom_vline(xintercept=attr(lib, "thresholds")[2], linetype="dotted", color="red") +
  geom_hline(yintercept=attr(genes, "thresholds")[1], linetype="dotted", color="blue") +
  geom_hline(yintercept=attr(genes, "thresholds")[2], linetype="dotted", color="blue")  
# flag the cells failing QC 
lungset$passQC <- ifelse(lungset$n.umi > attr(lib, "thresholds")[2] |
                           lungset$n.umi < attr(lib, "thresholds")[1] |
                           lungset$num_genes_expressed > attr(genes, "thresholds")[2] |
                           lungset$num_genes_expressed < attr(genes, "thresholds")[1] |
                           lungset$MT_perc > 25, #setting manually
                         FALSE, TRUE)
?table
table(lungset$passQC, lungset$clusters)

# For this particular dataset, the low quality cells all clustered together.
# Filtering based on hard or adaptive thresholds, or based on cluster membership
# both are valid options.

# Save!
?save_monocle_objects
save_monocle_objects(pbmc, "data/outputs/processed_pbmc_cds") #save for next class

################################################################################
###### Why are there are no rules or gold standards? ###########################
################################################################################
# Read in the metadata (column data) from 50K randomly subsampled barcodes from
# a benchmarking experiment our group did. These were performed with different
# technologies, on fixed samples, and much more shallowly sequenced.
df <- readRDS("data/inputs/benchmarking_metadata_50Kbarcodes_df.RDS")
str(df)

# Plot the distribution of UMIs and the threshold we determined above
# Does the red dotted line seem appropriate to you? HOPEFULLY NOT!
# Would you use the same threshold for all of these sample types? I SURE WOULD NOT!
ggplot(df, aes(x=paste(sample, assay), y=n.umi, fill=assay)) +
  facet_wrap(~samptype, scales="free_x", drop=T) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values=c("purple", "blue")) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=attr(lib, "thresholds")[1], linetype="dotted", color="red", size=1) #thresholds from the 10X PBMC data

# Plot the distribution of genes detected and the threshold we determined above
# Again, does this feel appropriate?
ggplot(df, aes(x=paste(sample, assay), y=num_genes_expressed, fill=assay)) +
  facet_wrap(~samptype, scales="free_x", drop=T) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values=c("purple", "blue")) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=attr(genes, "thresholds")[1], linetype="dotted", color="red", size=1) #thresholds from the 10X PBMC data

# Plot the distribution of MT and the threshold we determined above
# Again, does this feel appropriate?

# Because the dataset has both human and mouse derived samples, create a single
# column that we can plot on
df$MT_per <- ifelse(df$species == "hs", df$human_MT_reads/df$human_reads*100,
                    df$mouse_MT_reads/df$mouse_reads*100)

ggplot(df, aes(x=paste(sample, assay), y=MT_per, fill=assay)) +
  facet_wrap(~samptype, scales="free_x", drop=T) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  scale_fill_manual(values=c("purple", "blue")) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=attr(high.mito, "thresholds")[2], linetype="dotted", color="red", size=1) #thresholds from the 10X PBMC data

# Notice anything? Discuss.

################################################################################
###### Homework ################################################################
################################################################################      
# Play with thresholds
# Pick another dataset, download it, and play around with the above code. Or
# use your own data!
# Future tutorials will go more in depth, including multiple samples from 
# different runs, etc.

