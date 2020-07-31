library(fpop)


set.seed(1)
N <- 100
data.vec <- c(rnorm(N), rnorm(N, 2), rnorm(N))
fit <- Fpop(data.vec, N)
end.vec <- fit$t.est
change.vec <- end.vec[-length(end.vec)]
start.vec <- c(1, change.vec+1)
segs.list <- list()
for(seg.i in seq_along(start.vec)){
  start <- start.vec[seg.i]
  end <- end.vec[seg.i]
  seg.data <- data.vec[start:end]
  seg.mean <- mean(seg.data)
  segs.list[[seg.i]] <- data.frame(
    start, end,
    mean=seg.mean,
    seg.cost=sum((seg.data-seg.mean)^2))
}
segs <- do.call(rbind, segs.list)
plot(data.vec)
with(segs, segments(start-0.5, mean, end+0.5, mean, col="green"))
with(segs[-1,], abline(v=start-0.5, col="green", lty="dotted"))