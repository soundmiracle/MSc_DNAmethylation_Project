


###################################################################
#RESULT - CHART (For the ZxY table, show the number of rows (i.e number of 1-Mb of 100kb regions) where >Z dataset-pairs have r > Y (and r < -Y)
#Z={100,200,250,300}
#Y={0.2,0.5,0.8})

count.above.r=function(x,r) return(sum(x[x>r]))
count.below.r <- function(x, r) {
  return(sum(x <= r))
}

#1M cases
g=matrix(as.matrix(combined_matrix_1000000[,-c(1:4)]),ncol=dim(combined_matrix_1000000)[2]-4)
dim(g)
g[1,1:10]

apply(g,1,count.above.r,r=0.9)

summary(apply(g,1,count.above.r,r=0.9))

sum(apply(g,1,count.above.r,r=0.5)>200)
sum(apply(g,1,count.above.r,r=0.5)>250)
sum(apply(g,1,count.above.r,r=0.5)>300)
sum(apply(g,1,count.above.r,r=0.5)>350)

sum(apply(g,1,count.above.r,r=0.8)>200)
sum(apply(g,1,count.above.r,r=0.8)>250)
sum(apply(g,1,count.above.r,r=0.8)>300)
sum(apply(g,1,count.above.r,r=0.8)>350)

sum(apply(g,1,count.above.r,r=0.9)>200)
sum(apply(g,1,count.above.r,r=0.9)>250)
sum(apply(g,1,count.above.r,r=0.9)>300)
sum(apply(g,1,count.above.r,r=0.9)>350)

sum(apply(g,1,count.above.r,r=0.95)>200)

sum(apply(g,1,count.below.r,r=0)>200)
sum(apply(g,1,count.below.r,r=0)>250)

(1:dim(g)[1])[apply(g,1,count.above.r,r=0.9)>250]

#Find the area for making the Heatmap for data
combined_matrix_1000000[apply(g,1,count.above.r,r=0.9)>250,1:5]


combined_matrix_1000000[apply(g,1,count.below.r,r=0)>250,1:5]


sum(apply(g,1,count.below.r,r=0)>200)



#100K cases
f=matrix(as.matrix(combined_matrix_100k[,-c(1:4)]),ncol=dim(combined_matrix_100k)[2]-4)
dim(f)
f[1,1:10]

apply(f,1,count.above.r,r=0.9)

summary(apply(f,1,count.above.r,r=0.9))

sum(apply(f,1,count.above.r,r=0.5)>200)
sum(apply(f,1,count.above.r,r=0.5)>250)
sum(apply(f,1,count.above.r,r=0.5)>300)
sum(apply(f,1,count.above.r,r=0.5)>350)

sum(apply(f,1,count.above.r,r=0.8)>200)
sum(apply(f,1,count.above.r,r=0.8)>250)
sum(apply(f,1,count.above.r,r=0.8)>300)
sum(apply(f,1,count.above.r,r=0.8)>350)

sum(apply(f,1,count.above.r,r=0.9)>200)
sum(apply(f,1,count.above.r,r=0.9)>250)
sum(apply(f,1,count.above.r,r=0.9)>300)
sum(apply(f,1,count.above.r,r=0.9)>350)

sum(apply(f,1,count.above.r,r=0.95)>200)
sum(apply(f,1,count.above.r,r=0.95)>250)
sum(apply(f,1,count.above.r,r=0.95)>300)
sum(apply(f,1,count.above.r,r=0.95)>350)


sum(apply(f,1,count.below.r,r=0)>200)
sum(apply(f,1,count.below.r,r=0)>250)


pos_100k_0.95_250=combined_matrix_100k[apply(f,1,count.above.r,r=0.95)>250,1:4]

