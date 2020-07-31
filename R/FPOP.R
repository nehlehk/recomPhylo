library(fpop)
library(dplyr)
library(tidyr)

N <- 400
data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/500000/partial19.txt", character(), quote = "")
#data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/RAxML_perSiteLLs.likelihood_GTR", character(), quote = "")
data = noquote(data)
data = replace_na(data,0)
data = as.numeric(data)
data = data[4:500000]
#plot(data,main = "recombination rate= .00003 10000 ,branch length scaling factor = 0.2 and pen = 100")
plot(data,main = "node 19  - pen =400")
fit <- Fpop(data, N)
end.vec <- fit$t.est
print(end.vec)
abline(v= end.vec, col="blue")
  
  
  
  
  # print(mean(data[1:54176]))
  # print(var(data[1:54176]))
  # print(mean(data[54176:56997]))
  # print(var(data[54176:56997]))
  # print(mean(data[56997:71027]))
  # print(var(data[56997:71027]))
  # print(mean(data[71027:83301]))
  # print(var(data[71027:83301]))
  # print(mean(data[83301:499990]))
  # print(var(data[83301:499990]))
  
  
  # print(mean(data[1:22000]))
  # print(var(data[1:22000]))
  # print(mean(data[22000:25000]))
  # print(var(data[22000:25000]))
  # print(mean(data[25000:47000]))
  # print(var(data[25000:47000]))
  # print(mean(data[47000:50000]))
  # print(var(data[47000:50000]))
  # print(mean(data[50000:72000]))
  # print(var(data[50000:72000]))
  # print(mean(data[72000:75000]))
  # print(var(data[72000:75000]))
  # print(mean(data[75000:97000]))
  # print(var(data[75000:97000]))
  # print(mean(data[97000:100000]))
  # print(var(data[97000:100000]))
  
  
