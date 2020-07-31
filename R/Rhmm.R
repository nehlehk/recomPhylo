library('RcppHMM')

N = c("Low","Normal", "High")
A <- matrix(c(0.5, 0.3,0.2,
              0.2, 0.6, 0.2,
              0.1, 0.3, 0.6),  ncol= length(N), byrow = TRUE)

Mu <- matrix(c(0, 50, 100), ncol = length(N))
Sigma <- array(c(144, 400, 100), dim = c(1,1,length(N)))
Pi <- rep(1/length(N), length(N))
HMM.cont.univariate <- verifyModel(list( "Model"="GHMM",
                                         "StateNames" = N,
                                         "A" = A,
                                         "Mu" = Mu,
                                         "Sigma" = Sigma,
                                         "Pi" = Pi))


# Data simulation
set.seed(100)
length <- 100
observationSequence <- generateObservations(HMM.cont.univariate, length)

#Sequence decoding
hiddenStates <- viterbi(HMM.cont.univariate, observationSequence$Y)
print(hiddenStates)

