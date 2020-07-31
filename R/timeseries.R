# install.packages("https://CRAN.R-project.org/package=TTR")
# install.packages("TTR")


library(bcp)
library(dplyr)
library(tidyr)
library("TTR")

data = scan("/media/nehleh/295eaca0-f110-4b79-8cbe-bc99f9f61cbd/nehleh/0_Research/PhD/Data/simulationdata/recombination/clonalframe/RAxML_perSiteLLs.likelihood_GTR", character(), quote = "")
data = noquote(data)
data = replace_na(data,0)
data = data[1:100000]
x = as.numeric(data)

partialtimeseries <- ts(x)
SMA3 <- SMA(partialtimeseries, n= 10000)
E <- EMA(partialtimeseries, n= 3000 )


plot.ts(SMA3 , main="node12---n = 5000")
  
plot.ts(E ,   main="node12---n = 5000")

# a = acf(partialtimeseries, lag.max=100000 , plot=FALSE)
# pacf(partialtimeseries, lag.max=100000 )
