library(changepoint)

# Example of using the CROPS penalty in data set above
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))

out=cpt.mean(x, pen.value = c(4,1500),penalty = "CROPS",method = "PELT")

cpts.full(out)  # returns 7 segmentations for penalty values between 4 and 1500.  
# We find segmentations with 7, 5, 4, 3, 2, 1 and 0 changepoints. 
# Note that the empty final row indicates no changepoints.
pen.value.full(out) # gives associated penalty transition points
# CROPS does not give an optimal set of changepoints thus we may wish to explore further
plot(out,diagnostic=TRUE) 
# looks like the segmentation with 3 changepoints, 50,100,150 is the most appropriate
plot(out,ncpts=3) 
