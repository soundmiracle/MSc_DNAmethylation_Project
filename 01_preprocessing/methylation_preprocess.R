##### Methylation Preprocess


temp=commandArgs()
rds_files=as.character(temp[2])

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



