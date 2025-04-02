##### Local Methylation Correlation (LMC) Calculation


# Define the directory containing the RDS files
directory_path <- "/home/Datasets/"
save_path <- "/home/Corr_Result/"

# List all RDS files in the directory
##rds_files <- list.files(directory_path, pattern = "\\.RDS$", full.names = TRUE)

# Loop through each RDS file
for (rds_file in rds_files) {
  # Read the RDS file into the R environment
  cpg.mat <- readRDS(rds_file)
  # read in the data and store it into cpg.mat
  cpg.labels=rownames(cpg.mat)
  num.cpgs = dim(cpg.mat)[1]
  cpg.chromo=cpg.chromo.key[match(cpg.labels,cpg.key)]
  cpg.position=cpg.position.key[match(cpg.labels,cpg.key)]
  # write the code for sorting the cpg.chrome and cpg.position in order of numbers.
  cpg.chromo.position=cbind(cpg.labels, cpg.chromo, cpg.position)
  
  
  
  cpg.chromo.position[, 2] <- as.numeric(cpg.chromo.position[, 2])
  cpg.chromo.position[, 3] <- as.numeric(cpg.chromo.position[, 3])
  ordered_cpg.chromo.position <- cpg.chromo.position [order(as.numeric(cpg.chromo.position[, 2]), as.numeric(cpg.chromo.position[, 3])), ]
  ordered_cpg.label <- (ordered_cpg.chromo.position[, 1])
  ordered_cpg.chromo <- as.numeric(ordered_cpg.chromo.position[, 2])
  ordered_cpg.position <- as.numeric(ordered_cpg.chromo.position[, 3])
  ## find positions for each cpg in cpg.labels  --> cpg.position, cpg.chromo
  ## !!! might also need to order the rows of cpg.mat to be by chromosome and by increasing position !!!!
  
  
  
  
  order_index <- match(ordered_cpg.label, rownames(cpg.mat))
  ordered_cpg.mat <- cpg.mat[order_index, ]
  ordered_cpg.mat <- matrix(as.matrix(ordered_cpg.mat[,]),ncol=dim(ordered_cpg.mat)[2])
  
  
  
  ## calculate the correlation value btw cpgs from individuals with 10k range 
  
  num.comparisons=0
  for (i in 1:num.cpgs)
  {
    cpgs.tocompare=(1:num.cpgs)[ordered_cpg.chromo==ordered_cpg.chromo[i] & ordered_cpg.position> ordered_cpg.position[i] & ordered_cpg.position < (ordered_cpg.position[i]+10000)]
    num.comparisons=num.comparisons+length(cpgs.tocompare)
  }
  correlation.vec=distance.vec=cpg.label.1=cpg.label.2=rep(NA,num.comparisons)
  count=0
  for (i in 1:num.cpgs)
  {
    cpgs.tocompare=(1:num.cpgs)[ordered_cpg.chromo==ordered_cpg.chromo[i] & ordered_cpg.position> ordered_cpg.position[i] & ordered_cpg.position < (ordered_cpg.position[i]+10000)]
    if (length(cpgs.tocompare)>0)
    {
      for (j in 1:length(cpgs.tocompare))
      {
        distance.vec[(count+j)]=abs(ordered_cpg.position[i]-ordered_cpg.position[cpgs.tocompare[j]])
        correlation.vec[(count+j)]=round(cor(ordered_cpg.mat[i,],ordered_cpg.mat[cpgs.tocompare[j],], use = "complete.obs"),5)
        cpg.label.1[(count+j)]=ordered_cpg.label[i]
        cpg.label.2[(count+j)]=ordered_cpg.label[cpgs.tocompare[j]]
      }
    }
    count=count+length(cpgs.tocompare)
  }
  corr.result=cbind(cpg.label.1, cpg.label.2, correlation.vec, distance.vec)
  
  
  
  ## make the plot of r.sqr in 100 distance bins 
  
  #dist.bins=seq(0,10000,by=100)
  #mean.r.sqr.vec=rep(NA,length(dist.bins)-1)
  #for (i in 1:(length(dist.bins)-1))
  #{
  #  mean.r.sqr.vec[i]=mean(r.sqr.vec[distance.vec>dist.bins[i] & distance.vec<dist.bins[i+1]],na.rm=TRUE)
  #}
  
  # Construct the corresponding .RData file name
  rdata_file <- sub("\\.RDS$", ".RDS.RData", basename(rds_file))
  rdata_file_path <- file.path(save_path, rdata_file)
  
  # Save the environment to the .RData file ()
  save(corr.result, file = rdata_file_path) 
}