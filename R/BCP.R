library(bcp)
library(dplyr)
library(tidyr)
  
  
  
data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/500000-1/RAxML_perSiteLLs.likelihood_GTR", character(), quote = "")
data = noquote(data)
data = replace_na(data,0)
data = as.numeric(data)
data = data[4:500000]
bcp.0 <- bcp(data ,mcmc = 2000,burnin = 1)
plot(bcp.0, main="node***2000****12")

legacyplot(bcp.0)

#===================================================================================



x <- 1:100
b <- rep(c(3,-3), each=50)
y <- b*x + rnorm(100, sd=50)
plot(y)
bcp.3b <- bcp(y, x)
plot(bcp.3b, main="Linear Regression Change Point Example")
#===================================================================================


# 1 true change point at location 50; the predicting variable x is not related to location
x <- rnorm(100)
b <- rep(c(3,-3), each=50)
y <- b*x + rnorm(100)
plot(y)
bcp.3a <- bcp(y, x)


# in the two plots that follow, the location IDs are used as the plot characters
par(mfrow=c(1,2))
plot(y ~ x, type="n", main="Linear Regression: Raw Data")
text(x, y, as.character(1:100), col=(b/3)+2)
plot(y ~ x, type="n", main="Linear Regression: Posterior Means")
text(x, bcp.3a$posterior.mean[,1], as.character(1:100), col=(b/3)+2)
plot(bcp.3a, main="Linear Regression Change Point Example")
plot(bcp.3a, main="Linear Regression Change Point Example")
#===================================================================================
set.seed(5)
x <- rnorm(6, sd=3)
y <- rbind(cbind(rnorm(50, x[1]), rnorm(50, x[2]), rnorm(50, x[3])),
           cbind(rnorm(50, x[4]), rnorm(50, x[5]), rnorm(50, x[6])))

plot(y)
bcp.2a <- bcp(y)
plot(bcp.2a, main="Multivariate (k=3) Change Point Example")
plot(bcp.2a, separated=TRUE, main="Multivariate (k=3) Change Point Example")

#===================================================================================
set.seed(5)
x <- rep(c(0,1), each=50)
y <- x + rnorm(50, sd=1)
plot(y)
bcp.1b <- bcp(y)
plot(bcp.1b, main="Univariate Change Point Example")
legacyplot(bcp.1b)
#===================================================================================
set.seed(5)
x <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
plot(x)
bcp.1a <- bcp(x)
plot(bcp.1a, main="Univariate Change Point Example")
legacyplot(bcp.1a)

#===================================================================================

