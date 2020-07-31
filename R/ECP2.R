library("Rcpp")
library("ecp")
library(dplyr)
library(tidyr)

data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/12.txt", character(), quote = "")
data = noquote(data)
data = replace_na(data,0)
data = as.numeric(data)
data = data[20000:21000]

z = matrix(data, nrow = 1001 , ncol = 1 )
y= e.cp3o(z, K=10, minsize=100, alpha=1, verbose=FALSE)