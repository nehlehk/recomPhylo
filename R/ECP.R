library("Rcpp")
library("ecp")
library(dplyr)
library(tidyr)

data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/12.txt", character(), quote = "")
data = noquote(data)
data = replace_na(data,0)
data = as.numeric(data)
data = data[4:100000]

z = matrix(data, nrow = 99997 , ncol = 1 )
y= e.cp3o_delta(z, K=10, minsize=500, alpha=1, verbose=FALSE)


    