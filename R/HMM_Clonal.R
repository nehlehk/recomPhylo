library('RcppHMM')
library(dplyr)
library(tidyr)


data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/12.txt", character(), quote = "")
data = noquote(data)
data = replace_na(data,0)
data = as.numeric(data)
n = length(data)
avg = mean(data)
std = sd(data)


#states
N = c("Clonal", "non_Clonal")

#transProbs
A <- matrix(c(0.9999, 0.0001,
                0.00025, 0.99975),  ncol= length(N), byrow = TRUE)

#emissionProbs
Mu <- matrix(c(rnorm(n = 1 , mean = avg , sd = std), rnorm(n = 1 , mean = avg , sd = std)), ncol = length(N))

#covariance matrix of each state
Sigma <- array(c(1/std, 1/std), dim = c(1,1,length(N)))

#startProbs
Pi <- c(0.99, 0.01)

HMM.cont.univariate <- verifyModel(list( "Model"="GHMM",
                                         "StateNames" = N,
                                         "A" = A,
                                         "Mu" = Mu,
                                         "Sigma" = Sigma,
                                         "Pi" = Pi))




observationSequence <- matrix(c(data), ncol= n , byrow = TRUE )


newModel <- learnEM(HMM.cont.univariate,
                    observationSequence,
                    iter= 50,
                    delta = 1E-5,
                    print = FALSE)
print(newModel)


#Sequence decoding
hiddenStates <- viterbi(HMM.cont.univariate, observationSequence)
#print(hiddenStates)

fb <- forwardBackward(HMM.cont.univariate, observationSequence)


b = evaluation(HMM.cont.univariate, observationSequence, "b")
f = evaluation(HMM.cont.univariate, observationSequence, "f")

pos <- 1:100000
  
for (i in 1:n) {
  if (fb[i] =="Clonal") 
    y[i] = 2 
    else y[i] =1
}

plot(pos,y)


