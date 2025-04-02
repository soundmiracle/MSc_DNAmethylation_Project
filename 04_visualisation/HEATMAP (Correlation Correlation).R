##### RESULT - HEATMAP (Correlation Correlation)

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


load("C:/Users/mark9/Desktop/result/overall_result/Cor_Cor_Result1.RData")

# Load necessary libraries
library(ggplot2)
library(reshape2)



# Convert matrix row and column indices to a data frame with explicit numbering
cor_df <- melt(cor.mat)

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
  labs(title = "Heatmap of Correlation Matrix",
       x = "Dataset",
       y = "Dataset",
       fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(vjust = 0.5, hjust=1))

# Print the heatmap
print(heatmap_plot)