# Set the working directory to the folder where PLINK and input files are located
folder_path <- paste0("C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Colleagues/Dennis Lozada/plink_win64_20241022")
setwd(folder_path)

if (!require("tidyverse")) {install.packages("tidyverse", dependencies = TRUE)}
library(tidyverse)

# Define the PLINK executable
plink_executable <- "./plink.exe"  

# Run a basic PLINK command to check its version
plink_test_command <- paste(plink_executable, "--version")

# Execute the command
system(plink_test_command)

# Convert VCF SNP data to PLINK PED/MAP format
system("plink --vcf QV00490_112224_NMSU_Pepper_Pilot1_200X_allsamples_updated.Converted_unique.sorted.vcf --recode --out GBS192.plk")

#Convert PED/MAP to Binary Format
system("plink --ped GBS192.plk.ped --map GBS192.plk.map --make-bed --out GBS192.plk")

#LD-based SNP Pruning

# Pruning is important to reduce linkage disequilibrium and avoid bias in PCA.
# If your dataset has already been filtered to a small number of
# independent SNPs, you may skip this.
# Parameters:
# - 50 SNP window
# - 5 SNP step
# - r² threshold = 0.2

system("plink --bfile GBS192.plk --make-founders --indep-pairwise 50 5 0.2 --out pruned")

#Generate Pruned Dataset
system("plink --bfile GBS192.plk --extract pruned.prune.in --make-bed --out GBS192_pruned")

# Export Pruned Data as VCF (Optional)
system("plink --bfile GBS192_pruned --recode vcf --out GBS192_pruned_vcf")

#Calculate Genetic Distance Matrix
system("plink --allow-no-sex --nonfounders --bfile GBS192_pruned --distance-matrix --out dataForPCA")

# Load and Process Distance Matrix

dist_populations <- read.table("dataForPCA.mdist", header = FALSE)
fam <- data.frame(famids = read.table("dataForPCA.mdist.id")[, 1])
famInd <- data.frame(IID = read.table("dataForPCA.mdist.id")[, 2])

# Perform classical multidimensional scaling (PCA)
mds_populations <- cmdscale(dist_populations, eig = TRUE, k = 5)
eigenvec_populations <- cbind(fam, famInd, mds_populations$points)

# Calculate % variance explained by each PC
eigen_percent <- round((mds_populations$eig /
                          sum(mds_populations$eig)) * 100, 2)
print(eigen_percent)  # PC1, PC2, etc.

# Rename PCA columns for clarity
colnames(eigenvec_populations)[3:7] <- paste0("PC", 1:5)


#PCA plot without kmeans clustering

colnames(eigenvec_populations)[3:7] <- c("PC1", "PC2", "PC3", "PC4", "PC5")

# PCA Plot with samples colored by family ID (or any existing grouping)
ggplot(data = eigenvec_populations) +
  geom_point(aes(x = PC1, y = PC2, color = as.factor(famids)),
             size = 3, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(
    title = "PCA without K-means Clustering",
    x = paste0("PC1 (", eigen_percent[1], "% Explained)"),
    y = paste0("PC2 (", eigen_percent[2], "% Explained)")
  ) +
  theme_minimal(base_size = 14)

ggsave("PCA_no_clustering.png", width = 8, height = 6, dpi = 600, bg = "white")


# PCA Plot With K-means Clustering

# Perform k-means clustering using PC1 and PC2
set.seed(123)  # Reproducible results
k_clusters <- 7  # Define number of genetic groups expected
kmeans_result <- kmeans(eigenvec_populations[, c("PC1", "PC2")],
                        centers = k_clusters)

# Assign cluster labels
eigenvec_populations$Group <- paste0("Group ", kmeans_result$cluster)
eigenvec_populations$Group <- factor(eigenvec_populations$Group,
                                     levels = paste0("Group ", 1:k_clusters))

# Plot PCA with k-means cluster colors
ggplot(eigenvec_populations) +
  geom_point(aes(x = PC1, y = PC2, color = Group), size = 3) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(
    title = "PCA with K-means Clustering",
    x = paste0("PC1 (", eigen_percent[1], "% Explained)"),
    y = paste0("PC2 (", eigen_percent[2], "% Explained)"),
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14)

ggsave("PCA_kmeans_clustered.png", width = 8,
       height = 6, dpi = 600, bg = "white")


### NOTES

## PCA Visualization approaches in Genomic Studies
# PCA without Clustering
#-Shows natural genetic structure
#-Each point represents a genotype
#-Coloring by existing identifiers (e.g., family ID) reflects biological grouping but may show **overlap** due to admixture or shared ancestry.

## PCA with K-means Clustering
#-Assigns groups **objectively** based on similarity in principal component space (*unsupervised learning*).
#-Helps in grouping samples with **similar genetic makeup**, even if they originate from different biological sources.
#-Especially useful in genomics for:
#-Identifying subpopulations
#-Correcting for population structure** in GWAS
#-Visualizing diversity patterns** across samples

## Recommendation
#-Use K-means PCA plot when the goal is to identify or correct for hidden population structure (e.g., in GWAS or genomic prediction).
#-Use non-clustered PCA to examine raw diversity or visualize known groupings based on experimental design or pedigree.


