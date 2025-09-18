# set working directory 
setwd("")

# load data frame abiotic_data_fieldtrip_2023.csv
abiotic_data <- read.csv("abiotic_data_combined.csv", header = TRUE, sep = ",", dec = ".")

# show abiotic data
abiotic_data

# do a PCA based on the data in columns 2 to 10
abiotic_data_pca <- prcomp(abiotic_data[,4:12], scale = TRUE)

abiotic_data_pca

# Add the PCA scores to the data frame
pca_scores <- as.data.frame(abiotic_data_pca$x)
pca_scores

# Write PCA scores to a CSV file
write.csv(pca_scores, file = "pca_scores.csv", row.names = TRUE)

#pca_result$rotation
abiotic_data_pca$rotation
write.csv(abiotic_data_pca$rotation, "pca_loadings.csv", row.names = TRUE)

# Calculate the proportion of variance explained by each PC
explained_variance <- (abiotic_data_pca$sdev^2) / sum(abiotic_data_pca$sdev^2)

# Convert to percentages
explained_variance_percent <- explained_variance * 100

# Display the percentage of variance explained for PC1 and PC2
explained_variance_percent[1:2]


library(ggplot2)

#plot the PCA
plot(abiotic_data_pca$x[,1], abiotic_data_pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of abiotic data", pch = 19, col = "blue")

# add labels to the points
text(abiotic_data_pca$x[,1], abiotic_data_pca$x[,2], labels = rownames(abiotic_data), cex = 0.7, pos = 3)

# increase the size of the points in the plot
plot(abiotic_data_pca$x[,1], abiotic_data_pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of abiotic data", pch = 19, col = "blue", cex = 3)

# make triangles instead of points for all sites where year=2024
plot(abiotic_data_pca$x[,1], abiotic_data_pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of abiotic data", pch = ifelse(abiotic_data$year == 2024, 2, 19), col = "blue", cex = 3)


# set y- and x-axis limits to the same range
plot(abiotic_data_pca$x[,1], abiotic_data_pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of abiotic data", pch = ifelse(abiotic_data$year == 2024, 2, 19), col = "blue", cex = 3, xlim = c(-3.5, 6), ylim = c(-3.5, 6), cex.axis = 2, cex.lab = 2)



library(svglite)

# save the plot as a svg file
svglite("PCA_abiotic_data_combined1.svg", width = 10, height = 10)
P_abiotic <- plot(abiotic_data_pca$x[,1], abiotic_data_pca$x[,2], xlab = "PC1 (40.53%)", ylab = "PC2 (28.37%)", main = "PCA of abiotic data", pch = ifelse(abiotic_data$year == 2024, 2, 19), col = "blue", cex = 3, xlim = c(-3.5, 6), ylim = c(-3.5, 6), cex.axis = 2, cex.lab = 2)

text(abiotic_data_pca$x[,1], abiotic_data_pca$x[,2], labels = rownames(abiotic_data), cex = 0.7, pos = 3)

dev.off()


### heatmap PCA loadings

# Load necessary libraries
library(reshape2)

# Load PCA loadings data
data_loadings <- read.csv("pca_loadings.csv", header = TRUE, sep = ",", dec = ".")

# Extract the desired columns for PC1 and PC2
data <- data_loadings[, c("PC1", "PC2")]

# Reshape the data into long format using melt
data_long <- melt(as.matrix(data))

# Define custom labels for the x-axis (Var1)
# Assuming data_loadings has row names 1 to 9 mapped to specific variable names
lookup_labels <- c("Temperature", "pH", "DO", "Conductivity", "Phosphate", 
                   "Hardness", "Fe", "Cu", "NO2")

# Plot heatmap using ggplot

svglite("PCA_loadings.svg", width = 9, height = 4.5)
ggplot(data_long, aes(factor(Var1), Var2,  fill = value)) +  # Treat Var1 as a factor
  geom_tile() +  # Create tiles for the heatmap
  scale_fill_gradient(low = "yellow", high = "blue") +  # Define color gradient
  labs(x = "Parameter", y = "Principal Component", fill = "Loadings") +  # Define axis and legend labels
  theme_minimal() +  # Use a minimal theme for a clean look
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1),  # Increase x-axis text size
    axis.text.y = element_text(size = 18),  # Increase y-axis text size
    axis.title.x = element_text(size = 20, face = "bold"),  # Increase x-axis label size
    axis.title.y = element_text(size = 20, face = "bold"),  # Increase y-axis label size
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.title = element_text(size = 20, face = "bold")) +  
  scale_x_discrete(labels = lookup_labels)  # Map the numeric values to the descriptive labels

dev.off()
