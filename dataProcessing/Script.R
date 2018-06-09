library(data.table)
mydata<-fread("corrected_data_genome_score_31may18.txt")

pca1 = prcomp(mydata, scale = TRUE)
loadingsSquared = pca1$rotation^2
write.csv(loadingsSquared, file = "loadingsSquaredScore.csv")


pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  screeplot(x)
  screeplot(x,type="l")
  par(mfrow=c(1,1))
}

svdResult <- svd(mydata)
write.csv(svdResult$d, file="diagonalScore.csv")
write.csv(svdResult$u, file="uScore.csv")
write.csv(svdResult$v, file="vScore.csv")