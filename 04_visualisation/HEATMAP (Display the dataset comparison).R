###################################################################
#RESULT -  HEATMAP (display the dataset comparison)


# 1Mb - r > 0.9 - num.datasets>250- start.pos
start.pos.1M_0.9_250=combined_matrix_1000000[apply(g,1,count.above.r,r=0.9)>250,2]
# 1Mb - r > 0.9 - end.pos
end.pos.1M_0.9_250=combined_matrix_1000000[apply(g,1,count.above.r,r=0.9)>250,3]

unique.chrom=combined_matrix_1000000[apply(g,1,count.above.r,r=0.9)>250,1]
positions.chr=cpg.position.key[cpg.chromo.key==unique.chrom]
cpg.index.pos.1M_0.9_250=positions.chr[positions.chr>=start.pos.1M_0.9_250 & positions.chr<=end.pos.1M_0.9_250]

cpg.index.pos.1M_0.9_250=sort(cpg.index.pos.1M_0.9_250)
print(cpg.index.pos.1M_0.9_250)

cpg.index.1M_0.9_250=rep(NA,length(cpg.index.pos.1M_0.9_250))
for (b in 1:length(cpg.index.pos.1M_0.9_250)){
  cpg.index.1M_0.9_250[b]=cpg.key[cpg.position.key==cpg.index.pos.1M_0.9_250[b] & cpg.chromo.key==unique.chrom]
}
print(cpg.index.1M_0.9_250)






# Check each matrix in datasets
library(ggplot2)
library(dplyr)
library(gridExtra)
library(gplots)
library(dplyr)

results <- list()

for (i in 1:length(datasets)) {
  # Convert matrix to data frame
  df <- as.data.frame(datasets[[i]], stringsAsFactors = FALSE)
  colnames(df) <- c("cpg.1", "cpg.2", "r2")
  
  # Check if any CpG ID in cpg.index.1M_0.9_250 is in either cpg.1 or cpg.2 columns
  contains_cpg <- any(df$cpg.1 %in% cpg.index.1M_0.9_250) || any(df$cpg.2 %in% cpg.index.1M_0.9_250)
  
  # If it contains any CpG ID, add the index to the results
  if (contains_cpg) {
    results[[length(results) + 1]] <- list(matrix_index = i, data = df)
  }
}

# Print the results
if (length(results) > 0) {
  for (result in results) {
    cat("Matrix Index:", result$matrix_index, "\n")
    print(result$data)
  }
} else {
  cat("No matrices contain the specified CpG IDs.\n")
}



# Function to create a heatmap from a data frame and save it
create_and_save_heatmap <- function(df, title, filename) {
  # Prepare the data
  cpg_count <- length(cpg.index.1M_0.9_250)
  heatmap_matrix <- matrix(NA, nrow = cpg_count, ncol = cpg_count)
  rownames(heatmap_matrix) <- cpg.index.1M_0.9_250
  colnames(heatmap_matrix) <- cpg.index.1M_0.9_250
  
  for (i in 1:nrow(df)) {
    cpg1 <- df$cpg_label_1[i]
    cpg2 <- df$cpg_label_2[i]
    r2_value <- df$r_square_data[i]
    
    heatmap_matrix[cpg1, cpg2] <- r2_value
    heatmap_matrix[cpg2, cpg1] <- r2_value
  }
  
  # Create the heatmap plot
  heatmap_df <- as.data.frame(as.table(heatmap_matrix))
  colnames(heatmap_df) <- c("cpg1", "cpg2", "r2")
  
  p <- ggplot(heatmap_df, aes(x = cpg1, y = cpg2, fill = r2)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c("white", "red"),
                         na.value = "grey",
                         limits = c(0, 1),
                         breaks = c(0, 1),
                         labels = c("0", "1")) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +  # Remove the legend
    coord_fixed(ratio = 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title)
  
  # Save the heatmap
  ggsave(filename, plot = p, width = 5, height = 5)
  
  return(p)
}

# Create a list to store the heatmaps
heatmaps <- list()

# Generate and save heatmaps for datasets that contain CpG IDs in cpg.index.1M_0.9_250
for (i in 1:length(datasets)) {
  dataset_path <- names(datasets)[i]
  dataset_name <- sub("\\..*$", "", basename(dataset_path))  # Extract and clean the dataset name
  df <- as.data.frame(datasets[[i]], stringsAsFactors = FALSE)
  colnames(df) <- c("cpg_label_1", "cpg_label_2", "r_square_data")
  
  # Ensure the r2 column is numeric
  df$r_square_data <- as.numeric(df$r_square_data)
  
  # Filter the data
  filtered_data <- df %>%
    filter(cpg_label_1 %in% cpg.index.1M_0.9_250 & cpg_label_2 %in% cpg.index.1M_0.9_250)
  
  # If the filtered data is not empty, create and save a heatmap
  if (nrow(filtered_data) > 0) {
    heatmap_plot <- create_and_save_heatmap(filtered_data, paste(dataset_name), paste0(dataset_name, ".png"))
    heatmaps[[length(heatmaps) + 1]] <- heatmap_plot
  }
}

# Select random heatmaps to display in a 3x3 grid
set.seed(123)  # Set seed for reproducibility
selected_heatmaps <- sample(heatmaps, min(9, length(heatmaps)))

# Display selected heatmaps in a 3x3 grid
grid.arrange(grobs = selected_heatmaps, ncol = 3, nrow = 3)

###################################
#####100kb cases

# 100kb - r > 0.95 - num.datasets>250- start.pos
start.pos.100k_0.95_250=pos_100k_0.95_250[2,2]
# 100kb - r > 0.95 - num.datasets>250- end.pos
end.pos.100k_0.95_250=pos_100k_0.95_250[2,3]

unique.chrom=pos_100k_0.95_250[2,1]
positions.chr=cpg.position.key[cpg.chromo.key==unique.chrom]
cpg.index.pos.100k_0.95_250=positions.chr[positions.chr>=start.pos.100k_0.95_250 & positions.chr<=end.pos.100k_0.95_250]

cpg.index.pos.100k_0.95_250=sort(cpg.index.pos.100k_0.95_250)
print(cpg.index.pos.100k_0.95_250)

cpg.index.100k_0.95_250=rep(NA,length(cpg.index.pos.100k_0.95_250))
for (b in 1:length(cpg.index.pos.100k_0.95_250)){
  cpg.index.100k_0.95_250[b]=cpg.key[cpg.position.key==cpg.index.pos.100k_0.95_250[b] & cpg.chromo.key==unique.chrom]
}
print(cpg.index.100k_0.95_250)






# Check each matrix in datasets
library(ggplot2)
library(dplyr)
library(gridExtra)
library(gplots)
library(dplyr)

results <- list()

for (i in 1:length(datasets)) {
  # Convert matrix to data frame
  df <- as.data.frame(datasets[[i]], stringsAsFactors = FALSE)
  colnames(df) <- c("cpg.1", "cpg.2", "r2")
  
  # Check if any CpG ID in cpg.index.100k_0.95_250 is in either cpg.1 or cpg.2 columns
  contains_cpg <- any(df$cpg.1 %in% cpg.index.100k_0.95_250) || any(df$cpg.2 %in% cpg.index.100k_0.95_250)
  
  # If it contains any CpG ID, add the index to the results
  if (contains_cpg) {
    results[[length(results) + 1]] <- list(matrix_index = i, data = df)
  }
}

# Print the results
if (length(results) > 0) {
  for (result in results) {
    cat("Matrix Index:", result$matrix_index, "\n")
    print(result$data)
  }
} else {
  cat("No matrices contain the specified CpG IDs.\n")
}



# Function to create a heatmap from a data frame and save it
create_and_save_heatmap <- function(df, title, filename) {
  # Prepare the data
  cpg_count <- length(cpg.index.100k_0.95_250)
  heatmap_matrix <- matrix(NA, nrow = cpg_count, ncol = cpg_count)
  rownames(heatmap_matrix) <- cpg.index.100k_0.95_250
  colnames(heatmap_matrix) <- cpg.index.100k_0.95_250
  
  for (i in 1:nrow(df)) {
    cpg1 <- df$cpg_label_1[i]
    cpg2 <- df$cpg_label_2[i]
    r2_value <- df$r_square_data[i]
    
    heatmap_matrix[cpg1, cpg2] <- r2_value
    heatmap_matrix[cpg2, cpg1] <- r2_value
  }
  
  # Create the heatmap plot
  heatmap_df <- as.data.frame(as.table(heatmap_matrix))
  colnames(heatmap_df) <- c("cpg1", "cpg2", "r2")
  
  p <- ggplot(heatmap_df, aes(x = cpg1, y = cpg2, fill = r2)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c("white", "red"),
                         na.value = "grey",
                         limits = c(0, 1),
                         breaks = c(0, 1),
                         labels = c("0", "1")) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +  # Remove the legend
    coord_fixed(ratio = 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title)
  
  # Save the heatmap
  ggsave(filename, plot = p, width = 5, height = 5)
  
  return(p)
}

# Create a list to store the heatmaps
heatmaps <- list()

# Generate and save heatmaps for datasets that contain CpG IDs in cpg.index.100k_0.95_250
for (i in 1:length(datasets)) {
  dataset_path <- names(datasets)[i]
  dataset_name <- sub("\\..*$", "", basename(dataset_path))  # Extract and clean the dataset name
  df <- as.data.frame(datasets[[i]], stringsAsFactors = FALSE)
  colnames(df) <- c("cpg_label_1", "cpg_label_2", "r_square_data")
  
  # Ensure the r2 column is numeric
  df$r_square_data <- as.numeric(df$r_square_data)
  
  # Filter the data
  filtered_data <- df %>%
    filter(cpg_label_1 %in% cpg.index.100k_0.95_250 & cpg_label_2 %in% cpg.index.100k_0.95_250)
  
  # If the filtered data is not empty, create and save a heatmap
  if (nrow(filtered_data) > 0) {
    heatmap_plot <- create_and_save_heatmap(filtered_data, paste(dataset_name), paste0(dataset_name, ".png"))
    heatmaps[[length(heatmaps) + 1]] <- heatmap_plot
  }
}

# Select random heatmaps to display in a 3x3 grid
set.seed(123)  # Set seed for reproducibility
selected_heatmaps <- sample(heatmaps, min(9, length(heatmaps)))

# Display selected heatmaps in a 3x3 grid
grid.arrange(grobs = selected_heatmaps, ncol = 3, nrow = 3)