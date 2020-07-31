install.packages("TSPred", repo = 'https://cloud.r-project.org/bin/linux/ubuntu xenial/')

#Load TSPred package
  library("TSPred")

data <- c(1,2,1,2,10,2,1,2,1,2,3,4,5,6,2,5)


SlidingWindow("mean", data, 3, 1)