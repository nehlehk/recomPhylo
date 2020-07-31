library(changepoint)
library(dplyr)
library(tidyr)

data = scan("/home/nehleh/0_Research/PhD/Data/Ecoli-10_Aligned/likelihood_persite_GTR", character(), quote = "")
# data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/likelihood_GTR", character(), quote = "")



data = noquote(data)
data = replace_na(data,0)

data = as.numeric(data)

data = data[4:9640608]

out = cpt.mean(data,penalty='Manual',pen.value=5000,method='PELT',class = s = FALSE)

# cpt.mean(data,penalty="Manual",pen.value=1000,method="SegNeigh",class=FALSE) 

plot(data) 
abline(v=out, col="blue")

blocks = read.csv("/home/nehleh/0_Research/PhD/Data/xmfa/xmfa-blocks.csv", sep = ',')

myblock = blocks[,2]


abline(v=myblock, col="red")

print(out)
print(length(out))
  