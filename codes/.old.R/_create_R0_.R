R0 <- sapply(1:200, function(i) rnorm(n=100, mean=2, sd=.5))
write.csv(R0, file='data/vn_R0.csv', row.names=F)