##### High Correlation Region (HCR) Calculation 



options(scipen=999)

temp=commandArgs()

unique.chrom=as.character(temp[2])
block.size=as.double(temp[3])

out.file=paste("/home/Cor_Cor_result/BlocksofCpGs.corr.",unique.chrom,".",block.size,".txt",sep='')
out.file2=paste("/home /Cor_Cor_result/BlocksofCpGs.count.",unique.chrom,".",block.size,".txt",sep='')


# Define the directory containing the RData files
directory_path <- "/home/Corr_Result/" #######<= CHANGE the directory if available

save_path <- "/home/Cor_Cor/"


# List all RDS files in the directory
rdata_files <- list.files(directory_path, pattern = "\\.RData$", full.names = TRUE)
data.labels=rep(NA,length(rdata_files))
for (i in 1:length(data.labels)) data.labels[i]=strsplit(strsplit(rdata_files[i],split="Corr_Result/",fixed=T)[[1]][2],split=".RDS",fixed=T)[[1]][1]


####################################################################

#open the cpg.indicator file and load the cpg data
readfile=gzfile("/home/humanmethylation450_15017482_v1-2.csv.gz",open='r')


line2=scan(readfile,nlines=1,what='char',quiet=TRUE)  
num.cpgs.key=line.count=0
while(!is.na(line2[1]))
{
  if (sum(grep("cg",line2[1]))==1) num.cpgs.key=num.cpgs.key+1
  #if (sum(grep("cg",line2[1]))==0) print(c(line.count,line2[1]))
  line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
  line.count=line.count+1
}
close(readfile)
cpg.key=cpg.chromo.key=cpg.position.key=rep(NA,num.cpgs.key)
readfile=gzfile("/home/humanmethylation450_15017482_v1-2.csv.gz",open='r')


line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
cpg.count=1
while(!is.na(line2[1]))
{
  if (sum(grep("cg",line2[1]))==1)
  {
    cpg.key[cpg.count]=as.character(line3[1])
    cpg.chromo.key[cpg.count]=as.character(line3[12])
    cpg.position.key[cpg.count]=as.integer(line3[13])
    cpg.count=cpg.count+1
  }
  line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
  line3=strsplit(line2[1],split=',',fixed=T)[[1]]
}
close(readfile)


# Initialize a new list for storing R square data

datasets <- list()

cpgs.h=cpg.key[cpg.chromo.key==unique.chrom]
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
  corr.result = corr.result[is.element(cpg_label_1,cpgs.h),]
  
  
  datasets[[rdata_file]] <- corr.result
  # Rename the object based on the file name (without extension)
  #assign(x = tools::file_path_sans_ext(basename(rdata_file)), value = get(object_name, envir = env), envir = .GlobalEnv)
  
  # Optionally, clean up by removing the temporary environment
  rm(env)
  
  # clean up by removing the extra junks
  rm(corr.result, cpg_label_1, cpg_label_2, cor_data, r_square_data)
  
  
}


####################################################################
# #  correlation result for overall datasets
# 
# num.datasets=length(datasets)
# 
# cor.mat=count.mat=matrix(NA,nrow=num.datasets,ncol=num.datasets)
# for (j in 1:(num.datasets-1))
# {
#   dataset1.r2= as.numeric(datasets[[j]][,3])  ## column 3 of dataset1
#   dataset1.cpg.pairs=paste(datasets[[j]][,1],"-",datasets[[j]][,2],sep='')
#   for (k in (j+1):num.datasets)
#   {
#     dataset2.r2= as.numeric(datasets[[k]][,3])     ## column 3 of dataset2
#     dataset2.cpg.pairs=paste(datasets[[k]][,1],"-",datasets[[k]][,2],sep='')
#     overlap.cpg.pairs=intersect(dataset1.cpg.pairs,dataset2.cpg.pairs)
#     if (length(overlap.cpg.pairs)>0)
#     {
#       dataset1.r2.overlap=dataset1.r2[match(overlap.cpg.pairs,dataset1.cpg.pairs)]      ## vector containing correlation of CpGs that cpg i was compared to in dataset j that overlap with those compared in dataset k
#       dataset2.r2.overlap=dataset2.r2[match(overlap.cpg.pairs,dataset2.cpg.pairs)]      ## vector containing correlation of CpGs that cpg i was compared to in dataset k that overlap with those compared in dataset j
#       cor.mat[j,k]=cor.mat[k,j]=round(cor(dataset1.r2.overlap,dataset2.r2.overlap),3)
#       count.mat[j,k]=count.mat[k,j]=length(dataset1.r2.overlap)
#     }
#   }
# }
# 
# #to.print.cor=NULL
# #for (i in 1:(num.datasets-1)) to.print.cor=c(to.print.cor,round(cor.mat[i,(i+1):num.datasets],3))
# diag(cor.mat)=1
# to.print.cor=c(t(cor.mat))
# ## print out to.print.cor (i.e. cor.mat as a vector)
# ## print out count.mat as a vector into a separate file

####################################################################


#  correlation result for each cpgs

# Initialize the number of CpG sites and datasets
num.cpgs <- num.cpgs.key  # Total number of CpG sites on the array
num.datasets <- length(datasets)  # Total number of datasets

# Function to filter correlation results for a specific CpG
filter_cpg_correlation <- function(datasets, j, specific_cpg) {
  
  # Vector containing correlation of all CpGs that specific CpG was compared to in dataset j
  correlations <- as.numeric(datasets[[j]][,3])
  # CpG pairs in the format "CpG1-CpG2"
  cpg_pairs <- paste(datasets[[j]][,1], "-", datasets[[j]][,2], sep='')
  
  # Filter the results to include only pairs containing the specific CpG
  #specific_indices <- which(datasets[[j]][,1] == specific_cpg | datasets[[j]][,2] == specific_cpg)
  specific_indices = (1:length(datasets[[j]][,1]))[is.element(datasets[[j]][,1],specific_cpg) | is.element(datasets[[j]][,2],specific_cpg)]
  #match.firstpair=match(specific_cpg,datasets[[j]][,1])
  #match.secondpair=match(specific_cpg,datasets[[j]][,2])
  #specific_indices=sort(c(match.firstpair[!is.na(match.firstpair)],match.secondpair[!is.na(match.secondpair)]))
  
  
  filtered_correlations <- correlations[specific_indices]
  filtered_cpg_pairs <- cpg_pairs[specific_indices]
  
  
  return(list(correlations = filtered_correlations, cpg_pairs = filtered_cpg_pairs))
}

#cpg.datasets.cor.mat = cpg.datasets.count.mat = matrix(NA, nrow = num.cpgs, ncol = choose(num.datasets, 2))





header.val=NULL
for (i in 1:(length(data.labels)-1)) header.val=c(header.val,paste(data.labels[i],"-",data.labels[(i+1):length(data.labels)],sep=''))

# Set up the initial values

#unique.chrom=unique(cpg.chromo.key)
write(c("chromo","start.pos","end.pos","num.cpgs",header.val),file=out.file,ncolumns=length(header.val)+4)
write(c("chromo","start.pos","end.pos","num.cpgs",header.val),file=out.file2,ncolumns=length(header.val)+4)

# Loop through blocked CpG site
for (h in 1:length(unique.chrom))
{
  positions.chr=cpg.position.key[cpg.chromo.key==unique.chrom[h]]
  num.blocks=round((max(positions.chr)-min(positions.chr))/block.size,0)
  block.size.new=(max(positions.chr)-min(positions.chr))/num.blocks
  for (a in 1:num.blocks)
  {
    print(c(a,num.blocks))
    start.pos=min(positions.chr)+(a-1)*block.size.new
    end.pos=min(positions.chr)+a*block.size.new
    cpg.index.pos=positions.chr[positions.chr>=start.pos & positions.chr<=end.pos]
    
    if (length(cpg.index.pos)>10)
    {
      cpg.index=rep(NA,length(cpg.index.pos))
      for (b in 1:length(cpg.index.pos)){
        cpg.index[b]=cpg.key[cpg.position.key==cpg.index.pos[b] & cpg.chromo.key==unique.chrom[h]]
      }
      
      # Initialize matrices to store correlation values and counts
      cor.mat <- matrix(NA, nrow = num.datasets, ncol = num.datasets)
      count.mat <- matrix(NA, nrow = num.datasets, ncol = num.datasets)
      
      # Loop through each pair of datasets
      for (j in 1:(num.datasets - 1))
      {
        # Check if dataset1.r2 is not empty
        # Replace X with actual data for dataset j
        filter_dataset1=filter_cpg_correlation(datasets, j, cpg.index)
        dataset1.r2 <- as.numeric(filter_dataset1$correlations) # Vector containing correlation of all CpGs that CpG i was compared to in dataset j
        dataset1.cpgs <- filter_dataset1$cpg_pairs # Vector containing all CpGs that CpG i was compared to in dataset j
        
        if (length(dataset1.r2) > 0) 
        {
          for (k in (j + 1):num.datasets) 
          {
            print(c(a,j,k,num.datasets))
            # Replace Y with actual data for dataset k
            filter_dataset2=filter_cpg_correlation(datasets, k, cpg.index)
            dataset2.r2 <- as.numeric(filter_dataset2$correlations)  # Vector containing correlation of all CpGs that CpG i was compared to in dataset k
            dataset2.cpgs <- filter_dataset2$cpg_pairs  # Vector containing all CpGs that CpG i was compared to in dataset k
            
            # Check if dataset2.r2 is not empty
            if (length(dataset2.r2) > 0) {
              # Find overlapping CpG sites between dataset j and dataset k
              overlap.cpgs <- intersect(dataset1.cpgs, dataset2.cpgs)
              print(length(overlap.cpgs))
              if (length(overlap.cpgs)>1)
              {
                # Get the correlation values for the overlapping CpG sites
                dataset1.r2.overlap <- dataset1.r2[match(overlap.cpgs, dataset1.cpgs)]
                dataset2.r2.overlap <- dataset2.r2[match(overlap.cpgs, dataset2.cpgs)]
                
                # Calculate the correlation of the overlapping correlations and store in cor.mat
                cor.mat[j, k] <- cor.mat[k, j] <- round(cor(dataset1.r2.overlap, dataset2.r2.overlap), 3)
                
                # Store the count of overlapping CpG sites in count.mat
                count.mat[j, k] <- count.mat[k, j] <- length(dataset1.r2.overlap)
              }
            }
          }
        }
      }
      
      # Set the diagonal of the correlation matrix to 1 (correlation with itself)
      # diag(cor.mat) <- 1
      
      
      # Extract the upper triangle of the matrix excluding the diagonal
      upper_triangle_indices <- which(upper.tri(cor.mat))
      # Save the values in a vector
      to.print.cor <- cor.mat[upper_triangle_indices]
      to.print.count <- count.mat[upper_triangle_indices]
      
      
      # Print or save the correlation matrix vector (to.print.cor)
      #cpg.datasets.cor.mat[i,] = to.print.cor
      #cpg.datasets.count.mat[i,] = to.print.count
      
      write(c(unique.chrom[h],start.pos,end.pos,length(cpg.index.pos),to.print.cor),file=out.file,ncolumns=length(to.print.cor)+4,append=T)
      write(c(unique.chrom[h],start.pos,end.pos,length(cpg.index.pos),to.print.count),file=out.file2,ncolumns=length(to.print.cor)+4,append=T)
      
      # Here you can add code to save or print to.print.cor as needed
      # print(to.print.cor)
      
      # Optional: Save count.mat as well if needed
      # Here you can add code to save or print count.mat as needed
      # print(count.mat)
    }
  }
}
