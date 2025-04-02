###################################################################
#SAVE RESULT DATASETS
Datasets_path <- "/home/Datasets/"
txt_path <- "/home/Cor_Cor_result/" #######<= CHANGE the directory if available

save_path <- "/home/RESULT/"


# List all RDS files in the directory
rdata_files <- list.files(Datasets_path, pattern = "\\.RData$", full.names = TRUE)

data.labels=rep(NA,length(rdata_files))
for (i in 1:length(data.labels)) data.labels[i]=strsplit(strsplit(rdata_files[i],split="Corr_Result/",fixed=T)[[1]][2],split=".RDS",fixed=T)[[1]][1]

header.val=NULL
for (i in 1:(length(data.labels)-1)) header.val=c(header.val,paste(data.labels[i],"-",data.labels[(i+1):length(data.labels)],sep=''))

column.name=c('chromo', 'start.pos', 'end.pos', 'num.pairs',header.val)



# List all txt files in the directory
txt_files <- list.files(txt_path, pattern = "\\.txt$", full.names = TRUE)

result.labels=rep(NA,length(txt_files))
# Loop over each file path in txt_files
for (i in 1:length(txt_files)) {
  # Extract the file name without path and extension
  file_name <- strsplit(strsplit(txt_files[i], split = "result/", fixed = TRUE)[[1]][2], split = ".txt", fixed = TRUE)[[1]][1]
  result.labels[i] <- file_name
}

# Make the listset from txt files

cor.datasets <- list()

for (txt_file in txt_files) {
  # Read the file
  txt_result <- read.table(txt_file, header = TRUE, sep = " ", col.names = column.name, fill = TRUE)
  
  # Extract the base name and remove the extension
  file_name <- basename(txt_file)
  file_name <- sub("\\.txt$", "", file_name)
  
  # Save the result in the list with the cleaned file name
  cor.datasets[[file_name]] <- txt_result
}

# Construct the corresponding .RData file name
save_rdata_file <- "Cor_Cor_results.RData"
save_rdata_file_path <- file.path(save_path, save_rdata_file)

# Save the environment to the .RData file ()
save(cor.datasets, file = save_rdata_file_path) 



###################################################################
#RESULT - COMBINED_MATRIX_1M (Correlation each cpg-datasets)

# Filter the list to include only matrices with "1000000" or "corr" in their names
filtered_datasets_1000000 <- cor.datasets[grep("1000000.*corr|corr.*1000000", names(cor.datasets))]
# Combine the matrices row-wise
combined_matrix_1000000 <- do.call(rbind, filtered_datasets_1000000)


# Get the names of the columns to analyze (excluding 'chromo', 'start.pos', 'end.pos', 'num.pairs')
columns_to_analyze <- colnames(combined_matrix_1000000)[!(colnames(combined_matrix_1000000) %in% c('chromo', 'start.pos', 'end.pos', 'num.pairs'))]

# Calculate the average correlation for each column
average_correlations <- sapply(combined_matrix_1000000[, columns_to_analyze], mean, na.rm = TRUE)

# Use complete.cases to identify rows without NA values
complete_rows <- complete.cases(combined_matrix_1000000)

# Subset the matrix to exclude rows with NA values
combined_matrix_1000000 <- combined_matrix_1000000[complete_rows, ]



# Determine the column with the highest average correlation
max_average_correlation_column <- names(average_correlations)[which.max(average_correlations)]
max_average_correlation_value <- max(average_correlations)




###################################################################
#RESULT - COMBINED_MATRIX_100k (Correlation each cpg-datasets)

# Filter datasets to include only those with "100000" but not "1000000" and also include "corr"
filtered_datasets_100k <- cor.datasets[grep("corr.*100000(?!0)", names(cor.datasets), perl = TRUE)]

# Combine the matrices row-wise
combined_matrix_100k <- do.call(rbind, filtered_datasets_100k)


# Get the names of the columns to analyze (excluding 'chromo', 'start.pos', 'end.pos', 'num.pairs')
columns_to_analyze <- colnames(combined_matrix_100k)[!(colnames(combined_matrix_100k) %in% c('chromo', 'start.pos', 'end.pos', 'num.pairs'))]

# Calculate the average correlation for each column
average_correlations_100k <- sapply(combined_matrix_100k[, columns_to_analyze], mean, na.rm = TRUE)

# Use complete.cases to identify rows without NA values
complete_rows <- complete.cases(combined_matrix_100k)

# Subset the matrix to exclude rows with NA values
combined_matrix_100k <- combined_matrix_100k[complete_rows, ]


###################################################################
#RESULT - HEATMAP (Correlation each cpg-datasets in 1Mb or 100kb gate)

# Filter the list to include only matrices with "1000000" or "corr" in their names
filtered_datasets_1000000 <- cor.datasets[grep("1000000.*corr|corr.*1000000", names(cor.datasets))]

# Combine the matrices row-wise
combined_matrix_1000000 <- do.call(rbind, filtered_datasets_1000000)
# Rounding the start.pos and end.pos
combined_matrix_1000000[,2]=round(combined_matrix_1000000[,2])
combined_matrix_1000000[,3]=round(combined_matrix_1000000[,3])

# Use complete.cases to identify rows without NA values
complete_rows <- complete.cases(combined_matrix_1000000)

# Subset the matrix to exclude rows with NA values
combined_matrix_1000000 <- combined_matrix_1000000[complete_rows, ]


aver_1000000=rep(NA, 351)
for (i in 351){
  aver_1000000_value=mean(as.numeric(combined_matrix_1000000[,i+4]))
  aver_1000000[i]=aver_1000000_value
}

# Number of combinations of 27 choose 2
num_combinations <- choose(27, 2)

# Initialize the vector to store average values
aver_1000000 <- rep(NA, num_combinations)


# Loop through each combination and calculate the mean of the corresponding column in the matrix
for (i in 1:num_combinations) {
  column_index <- i + 4  # Adjusting the index to start from the 5th column
  if (column_index <= ncol(combined_matrix_1000000)) {
    aver_1000000[i] <- mean(as.numeric(combined_matrix_1000000[, column_index]), na.rm = TRUE)
  }
}

# Initialize a 27x27 matrix with NA values
cor_mat_1000000 <- matrix(NA, nrow = 27, ncol = 27)

# Fill the upper triangular part of the matrix
index <- 1
for (i in 1:26) {
  for (j in (i + 1):27) {
    cor_mat_1000000[i, j] <- aver_1000000[index]
    index <- index + 1
  }
}

# Mirror the upper triangular part to the lower triangular part
for (i in 1:27) {
  for (j in 1:27) {
    if (is.na(cor_mat_1000000[j, i])) {
      cor_mat_1000000[j, i] <- cor_mat_1000000[i, j]
    }
  }
}

# Optional: Set the diagonal to 1 (as the correlation of each dataset with itself is 1)
diag(cor_mat_1000000) <- 1






# Load necessary libraries
library(ggplot2)
library(reshape2)



# Convert matrix row and column indices to a data frame with explicit numbering
cor_df <- melt(cor_mat_1000000)

# Dataset names without .RDS extension
dataset_names <- c("Blood_Cauc", "Blood_Hisp", "Blood_Japan", "Blood_Mexican", "Blood_PuertoRican", 
                   "Buccals_Cauc", "Buccals_Sing_9mo", "BulkFrontalCortex", "CD4+_Estonian", "CD8+_Estonian",
                   "CordBlood", "PBMC_AfricanAmerican", "Placenta", "Saliva_Cauc", "Saliva_Hisp", 
                   "Skin_UKTwin", "SuperiorTemporalGyrus", "TCGA_bladder", "TCGA_breast", "TCGA_colon", 
                   "TCGA_head_and_neck", "TCGA_kidney", "TCGA_liver", "TCGA_lung", "TCGA_prostate", 
                   "TCGA_thyroid", "TCGA_uterus")

# Add the dataset names as factors
cor_df$Var1 <- factor(cor_df$Var1, levels = 1:27, labels = dataset_names)
cor_df$Var2 <- factor(cor_df$Var2, levels = 1:27, labels = dataset_names)

# Create the heatmap
heatmap_plot <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Heatmap of Correlation Matrix (1Mb gate)",
       x = "Dataset",
       y = "Dataset",
       fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(vjust = 0.5, hjust=1))

# Print the heatmap
print(heatmap_plot)




## 100kb cases 

# Filter datasets to include only those with "100000" but not "1000000" and also include "corr"
filtered_datasets_100k <- cor.datasets[grep("corr.*100000(?!0)", names(cor.datasets), perl = TRUE)]

# Combine the matrices row-wise
combined_matrix_100k <- do.call(rbind, filtered_datasets_100k)


# Get the names of the columns to analyze (excluding 'chromo', 'start.pos', 'end.pos', 'num.pairs')
columns_to_analyze <- colnames(combined_matrix_100k)[!(colnames(combined_matrix_100k) %in% c('chromo', 'start.pos', 'end.pos', 'num.pairs'))]

# Calculate the average correlation for each column
average_correlations_100k <- sapply(combined_matrix_100k[, columns_to_analyze], mean, na.rm = TRUE)

# Use complete.cases to identify rows without NA values
complete_rows <- complete.cases(combined_matrix_100k)

# Subset the matrix to exclude rows with NA values
combined_matrix_100k <- combined_matrix_100k[complete_rows, ]


aver_100k=rep(NA, 351)
for (i in 351){
  aver_100k_value=mean(as.numeric(combined_matrix_100k[,i+4]))
  aver_100k[i]=aver_100k_value
}

# Number of combinations of 27 choose 2
num_combinations <- choose(27, 2)

# Initialize the vector to store average values
aver_100k <- rep(NA, num_combinations)


# Loop through each combination and calculate the mean of the corresponding column in the matrix
for (i in 1:num_combinations) {
  column_index <- i + 4  # Adjusting the index to start from the 5th column
  if (column_index <= ncol(combined_matrix_100k)) {
    aver_100k[i] <- mean(as.numeric(combined_matrix_100k[, column_index]), na.rm = TRUE)
  }
}

# Initialize a 27x27 matrix with NA values
cor_mat_100k <- matrix(NA, nrow = 27, ncol = 27)

# Fill the upper triangular part of the matrix
index <- 1
for (i in 1:26) {
  for (j in (i + 1):27) {
    cor_mat_100k[i, j] <- aver_100k[index]
    index <- index + 1
  }
}

# Mirror the upper triangular part to the lower triangular part
for (i in 1:27) {
  for (j in 1:27) {
    if (is.na(cor_mat_100k[j, i])) {
      cor_mat_100k[j, i] <- cor_mat_100k[i, j]
    }
  }
}

# Optional: Set the diagonal to 1 (as the correlation of each dataset with itself is 1)
diag(cor_mat_100k) <- 1






# Load necessary libraries
library(ggplot2)
library(reshape2)



# Convert matrix row and column indices to a data frame with explicit numbering
cor_df <- melt(cor_mat_100k)

# Dataset names without .RDS extension
dataset_names <- c("Blood_Cauc", "Blood_Hisp", "Blood_Japan", "Blood_Mexican", "Blood_PuertoRican", 
                   "Buccals_Cauc", "Buccals_Sing_9mo", "BulkFrontalCortex", "CD4+_Estonian", "CD8+_Estonian",
                   "CordBlood", "PBMC_AfricanAmerican", "Placenta", "Saliva_Cauc", "Saliva_Hisp", 
                   "Skin_UKTwin", "SuperiorTemporalGyrus", "TCGA_bladder", "TCGA_breast", "TCGA_colon", 
                   "TCGA_head_and_neck", "TCGA_kidney", "TCGA_liver", "TCGA_lung", "TCGA_prostate", 
                   "TCGA_thyroid", "TCGA_uterus")

# Add the dataset names as factors
cor_df$Var1 <- factor(cor_df$Var1, levels = 1:27, labels = dataset_names)
cor_df$Var2 <- factor(cor_df$Var2, levels = 1:27, labels = dataset_names)

# Create the heatmap
heatmap_plot <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Heatmap of Correlation Matrix (100kb gate)",
       x = "Dataset",
       y = "Dataset",
       fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(vjust = 0.5, hjust=1))

# Print the heatmap
print(heatmap_plot)
