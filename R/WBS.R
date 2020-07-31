library("wbs")

x = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/500000/RAxML_perSiteLLs.likelihood_GTR", character(), quote = "")

#x <- rnorm(300) + c(rep(1,50),rep(0,250))

#plot(x , main = "WBS")

w <- wbs(x)
w.cpt <- changepoints(w,th = 2 )
w.cpt
plot(w)



s <- sbs(x)
s.cpt <- changepoints(s , th = 3 )
s.cpt
plot(s)




