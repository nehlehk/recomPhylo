library(cpm)

x <- scan("/home/nehleh/0_Research/PhD/Data/simulationdata/likelihood_GTR", character(), quote = "")

result <- processStream(x,"Kolmogorov-Smirnov",startup=1000, ARL0 = 10000 ,lambda=0.3)

plot(x)

for (i in 1:length(result$changePoints)) {
abline(v=result$changePoints[i], lty=2 ,col="blue" ,lwd=1)
}