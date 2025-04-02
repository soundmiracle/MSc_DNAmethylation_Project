###################################################################
#RESULT -  violin plot (display the dataset comparison)


# Load necessary packages
library(ggplot2)
library(reshape2)
library(dplyr)

rds_file = "C:/Users/mark9/Desktop/Datasets/Saliva_Hisp.RDS"
cpg.mat <- readRDS(rds_file)

# Find the intersection of cpg.index.100k_0.95_250 vector and row names of cpg.mat, then sort rows accordingly
common_cpg_ids <- intersect(cpg.index.100k_0.95_250, rownames(cpg.mat))
filtered.cpg.mat <- cpg.mat[common_cpg_ids, ]yyouyout

# Sort the rows of filtered.cpg.mat in the order of cpg.index.100k_0.95_250
filtered.cpg.mat <- filtered.cpg.mat[match(cpg.index.100k_0.95_250, common_cpg_ids), ]

# Convert the beta scores to a data frame
filtered.cpg.df <- as.data.frame(filtered.cpg.mat)
filtered.cpg.df$CpG_ID <- factor(rownames(filtered.cpg.df), levels = cpg.index.100k_0.95_250)

# Convert the data frame to long format
melted.df <- melt(filtered.cpg.df, id.vars = "CpG_ID")

# Calculate the mean and standard deviation for each CpG (excluding NA values)
stats.df <- melted.df %>%
  group_by(CpG_ID) %>%
  summarize(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))

# Create the violin plot
ggplot(melted.df, aes(x = CpG_ID, y = value)) +
  geom_violin(trim = FALSE) +
  geom_point(data = stats.df, aes(x = CpG_ID, y = mean), color = "red") +
  geom_errorbar(data = stats.df, aes(x = CpG_ID, ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "blue") +
  labs(x = "CpG ID", y = "Beta Score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))