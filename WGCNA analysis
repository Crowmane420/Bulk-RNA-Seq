# Load required libraries
library(tidyverse)  # For data manipulation and visualization
library(magrittr)   # For pipe operators
library(WGCNA)      # For weighted gene co-expression network analysis
library(DESeq2)     # For differential expression analysis
allowWGCNAThreads() # Allow multithreading in WGCNA

# ==== Load and clean data
# Load gene count data with gene IDs as row names
# Ensure the file path is correct for your system
data_WGCNA <- read.csv("/home/goutham/Desktop/RNA_Seq/Phoenix/data/salmon.merged.gene_counts.csv", row.names = 1)

# Column selection for pivot_longer
col_sel <- names(data_WGCNA)

# Load metadata for sample groups
mdata <- read.csv("/home/goutham/Desktop/RNA_Seq/Phoenix/meta/meta.csv", row.names = 1)

# Create a mapping of sample names to groups
group_map <- mdata$Groups
names(group_map) <- rownames(mdata) 

# Transform data into long format and assign groups
mdata_long <- data_WGCNA %>%
  pivot_longer(
    cols = all_of(col_sel),
    names_to = "name",
    values_to = "value"
  ) %>%
  mutate(
    group = group_map[name]  # Add group information
  )

# Plot RNA-Seq counts grouped by sample groups
p <- mdata_long %>%
  ggplot(aes(x = name, y = value)) +
  geom_violin() +
  geom_point(alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Treatment Groups", y = "RNA Seq Counts") +
  facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")

# Print the plot
print(p)

# ==== Prepare data for DESeq2 analysis
# Convert gene count data to a matrix
de_input <- as.matrix(data_WGCNA)

# Create metadata for DESeq2
meta_df <- data.frame(Sample = colnames(de_input)) %>%
  mutate(Type = gsub("_REP.*", "", Sample))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(de_input),  # Round counts to integers
  colData = meta_df,           # Metadata
  design = ~ Type              # Experimental design formula
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Perform variance stabilizing transformation (VST)
vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- assay(vsd)  # Extract normalized data

# Calculate row variances
rv_wpn <- rowVars(wpn_vsd)

# Filter highly variable genes using the 95th percentile
q95_wpn <- quantile(rv_wpn, 0.95)
expr_normalized <- wpn_vsd[rv_wpn > q95_wpn, ]

# Transform normalized data to long format
expr_normalized_df <- as.data.frame(expr_normalized) %>%
  mutate(Gene_id = row.names(expr_normalized)) %>%
  pivot_longer(-Gene_id, names_to = "name", values_to = "value")

# ==== Visualize normalized expression data
expr_normalized_df %>%
  ggplot(aes(x = name, y = value)) +
  geom_violin() +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95th Quantile Expression",
    x = "Treatment",
    y = "Normalized Expression"
  )

# Add group information to normalized data
meta_df <- data.frame(Sample = colnames(de_input)) %>%
  mutate(Group = gsub("_REP.*", "", Sample))

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(Gene_id = row.names(expr_normalized)) %>%
  pivot_longer(-Gene_id, names_to = "name", values_to = "value") %>%
  left_join(meta_df, by = c("name" = "Sample"))

# Plot grouped normalized expression
grouped_plot <- expr_normalized_df %>%
  ggplot(aes(x = name, y = value, fill = Group)) +
  geom_violin() +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ Group, scales = "free_x") +
  ylim(0, NA) +
  labs(
    title = "Grouped Normalized Expression",
    x = "Replicates",
    y = "Normalized Expression"
  )

# Print the grouped plot
print(grouped_plot)

# ==== Perform WGCNA
# Transpose normalized data for WGCNA
input_mat <- t(expr_normalized)

# Determine soft threshold for network construction
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(input_mat, powerVector = powers, verbose = 5)

# Plot scale-free topology model fit and mean connectivity
par(mfrow = c(1, 2))
cex1 <- 0.9

# Plot scale independence
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale Independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

# Plot mean connectivity
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

# Select the power with the best fit
picked_power <- 9

# Construct network using WGCNA
cor <- WGCNA::cor  # Use WGCNA's correlation function
netwk <- blockwiseModules(input_mat,
                          power = picked_power,
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = FALSE,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "ER",
                          numericLabels = TRUE,
                          verbose = 3)

# Assign module colors
mergedColors <- labels2colors(netwk$colors)

# Plot dendrogram with module colors
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

# Export gene modules
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors))
write_delim(module_df, file = "gene_modules.txt", delim = "\t")

# Calculate module eigengenes
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)

# Prepare module-trait relationships for plotting
module_order <- names(MEs0) %>% gsub("ME", "", .)
MEs0$treatment <- row.names(MEs0)
mME <- MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order))

# Plot module-trait relationships
mME %>% ggplot(aes(x = treatment, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Module-Trait Relationships", y = "Modules", fill = "Correlation")

# ==== Generate edge list
# Select modules of interest
modules_of_interest <- c("grey", "turquoise", "brown", "yellow", "blue")
submod <- module_df %>% subset(colors %in% modules_of_interest)

# Extract expression data for genes in modules of interest
row.names(module_df) <- module_df$gene_id
subexpr <- expr_normalized[submod$gene_id, ]

# Prepare data for edge list
tom <- TOMsimilarityFromExpr(t(subexpr), power = picked_power)
row.names(tom) <- row.names(subexpr)
colnames(tom) <- row.names(subexpr)
edge_list <- data.frame(tom) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1) %>%
  rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1 == gene2)) %>%
  mutate(
    module1 = module_df[gene1, ]$colors,
    module2 = module_df[gene2, ]$colors)

# Export edge list
write_delim(edge_list, file = "edgelist.tsv", delim = "\t")
