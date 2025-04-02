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
#RESULT - HISTOGRAM (Correlation each cpg-datasets in 1000000 gate)

# Filter the list to include only matrices with "1000000" or "corr" in their names
filtered_datasets_1000000 <- cor.datasets[grep("1000000.*corr|corr.*1000000", names(cor.datasets))]

# Combine the matrices row-wise
combined_matrix_1000000 <- do.call(rbind, filtered_datasets_1000000)

# Get the names of the columns to analyze (excluding 'chromo', 'start.pos', 'end.pos', 'num.pairs')
columns_to_analyze <- colnames(combined_matrix_1000000)[!(colnames(combined_matrix_1000000) %in% c('chromo', 'start.pos', 'end.pos', 'num.pairs'))]


# Calculate the average r-squared for each column
average_correlations <- sapply(combined_matrix_1000000[, columns_to_analyze], mean, na.rm = TRUE)

# Determine the column with the highest average r-squared
max_average_correlation_column <- names(average_correlations)[which.max(average_correlations)]
max_average_correlation_value <- max(average_correlations)

# Determine the column with the lowest average r-squared
min_average_correlation_column <- names(average_correlations)[which.min(average_correlations)]
min_average_correlation_value <- min(average_correlations)

# Plot histogram for the column with the highest average r-squared
data_vector <- combined_matrix_1000000[[max_average_correlation_column]]
hist(data_vector, main = paste("Histogram of", max_average_correlation_column), xlab = "R-squared Values", col = "blue", border = "black")

# Plot histogram for the column with the lowest average r-squared
data_vector <- combined_matrix_1000000[[min_average_correlation_column]]
hist(data_vector, main = paste("Histogram of", min_average_correlation_column), xlab = "R-squared Values", col = "red", border = "black")



# Plot density for the column with the highest average r-squared
data_vector_high <- combined_matrix_1000000[[max_average_correlation_column]]
density_high <- density(data_vector_high, na.rm = TRUE)

# Plot density for the column with the lowest average r-squared
data_vector_low <- combined_matrix_1000000[[min_average_correlation_column]]
density_low <- density(data_vector_low, na.rm = TRUE)

# Determine the maximum y value to set the y-axis limit
max_y <- max(c(density_high$y, density_low$y))

# Set up the plot with adjusted x and y-axis limits
plot(density_high, main = "Density Plot of R-squared Values", xlab = "R-squared Values", ylab = "Density", col = "blue", lwd = 2, ylim = c(0, max_y * 1.1), xlim = c(0, 1))
lines(density_low, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c(max_average_correlation_column, min_average_correlation_column), col = c("blue", "red"), lwd = 2)





# Plot histogram for the average r-squared values
hist(average_correlations, main = "Histogram of Average R-squared Values", xlab = "Average R-squared Values", col = "green", border = "black")


# Set up the plot area
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 40)) # Adjust x and y limits as necessary
title(main = "Overlapping Density Plots", xlab = "Correlation Values", ylab = "Density")

# Define colors for the plots
colors <- rainbow(length(columns_to_analyze))

# Plot density for each column
for (i in seq_along(columns_to_analyze)) {
  col <- columns_to_analyze[i]
  data_vector <- combined_matrix_1000000[[col]]
  
  # Plot density
  lines(density(data_vector, na.rm = TRUE), col = colors[i])
}

