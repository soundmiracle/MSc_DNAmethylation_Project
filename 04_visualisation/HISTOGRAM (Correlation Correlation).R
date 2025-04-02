##### RESULT - HISTOGRAM (Correlation Correlation)


# Define the directory containing the RData files
directory_path <- "C:/Users/mark9/Desktop/result/" #######<= CHANGE the directory if available

save_path <- "C:/Users/mark9/Desktop/"

# List all RDS files in the directory
rdata_files <- list.files(directory_path, pattern = "\\.RData$", full.names = TRUE)


# Initialize a new list for storing R square data

datasets <- list()

for (rdata_file in rdata_files)
{
  # Create a new environment for loading the file
  env <- new.env()
  
  # Load the file into the new environment
  load(rdata_file, envir = env)
  
  
  # Get the name of the first object loaded (assuming one main object per file)
  object_name <- ls(env)[1]
  corr.result <- env$corr.result
  cor_data <-as.numeric(corr.result[,3])  # Assuming the third column is R square
  r_square_data <- cor_data^2
  cpg_label_1 = corr.result[,1]
  cpg_label_2 = corr.result[,2]
  corr.result = cbind(cpg_label_1, cpg_label_2, r_square_data)
  
  
  
  datasets[[rdata_file]] <- corr.result
  # Rename the object based on the file name (without extension)
  assign(x = tools::file_path_sans_ext(basename(rdata_file)), value = get(object_name, envir = env), envir = .GlobalEnv)
  
  # Optionally, clean up by removing the temporary environment
  rm(env)
  
  # clean up by removing the extra junks
  rm(corr.result, cpg_label_1, cpg_label_2, cor_data, r_square_data)
  
  
}



###################################################################


load(rdata_file, envir = env)


# Load necessary libraries
library(ggplot2)
library(dplyr)

# List of dataset names without .RDS extension
dataset_names <- c("Blood_Cauc", "Blood_Hisp", "Blood_Japan", "Blood_Mexican", "Blood_PuertoRican", 
                   "Buccals_Cauc", "Buccals_Sing_9mo", "BulkFrontalCortex", "CD4+_Estonian", "CD8+_Estonian",
                   "CordBlood", "PBMC_AfricanAmerican", "Placenta", "Saliva_Cauc", "Saliva_Hisp", 
                   "Skin_UKTwin", "SuperiorTemporalGyrus", "TCGA_bladder", "TCGA_breast", "TCGA_colon", 
                   "TCGA_head_and_neck", "TCGA_kidney", "TCGA_liver", "TCGA_lung", "TCGA_prostate", 
                   "TCGA_thyroid", "TCGA_uterus")

# Function to calculate mean r-squared values for each distance bin
calculate_mean_r_squared <- function(dataset_name) {
  # Assuming each dataset is loaded into a variable with its name suffixed by _RDS
  dataset <- get(paste0(dataset_name, ".RDS"))
  
  distance.vec <- as.numeric(dataset[,4])
  r.sqr.vec <- as.numeric(dataset[,3])
  dist.bins <- seq(0, 10000, by = 100)
  mean.r.sqr.vec <- rep(NA, length(dist.bins) - 1)
  
  for (i in 1:(length(dist.bins) - 1)) {
    mean.r.sqr.vec[i] <- mean(r.sqr.vec[distance.vec > dist.bins[i] & distance.vec < dist.bins[i + 1]], na.rm = TRUE)
  }
  
  data.frame(
    DistanceBin = seq(50, 9950, by = 100),
    AvgRSqrValue = mean.r.sqr.vec,
    Dataset = dataset_name
  )
}

# Combine results from all datasets into a single data frame
combined_results <- bind_rows(lapply(dataset_names, calculate_mean_r_squared))

# Generate the plot
ggplot(combined_results, aes(x = DistanceBin, y = AvgRSqrValue, color = Dataset)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    x = "Distance bins in 100 base pair increments",
    y = "Average r-squared value",
    color = "Dataset"
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.title = element_blank())
