

## Effect of hierarchical levels on variance of final level


HGammaVar <- function(K, shape, mean){
  
  AdjustShape <- function(shape, K){
    shape * K + ((K - 1)/(shape+1))
  }
  
  scale <- mean / shape
  
  shape <- AdjustShape(shape, K)
  
  K <- K+1
  
  g <- rep(NA, K)
  g[1] <- mean 
  
  for (i in 2:K){
    g[i] <- rgamma(1, shape = shape, scale = g[i-1] / shape) 
  }
  
  return(g[K])
  
}

ExpVar <- function(shape, mean){
  
  scale <- mean/(shape)
  var1 <- (shape * scale^2)
  var1 
}


shp <- 1.5


x1 <- replicate(100^2, HGammaVar(K = 8, shape = shp, mean = 1))
x1 <- matrix(x1, ncol = 100)
colMeans(x1)
x1.var <- apply(x1, 2, var)
mean(x1.var)
hist(x1.var, 25, xlim = c(0, max(x1.var)))

abline(v = ExpVar(shape = shp, mean = 1), col = "Red")
abline(v = mean(x1.var), col = "Blue")
abline(v = median(x1.var), col = "Green")




### With Gaussian

HGaus <- function(K, sd, mean){
  
  K <- K+1
  
  g <- rep(NA, K)
  g[1] <- mean 
  
  for (i in 2:K){
    g[i] <- rnorm(1, mean = g[i-1], sd = sd) 
  }
  
  return(g[K])
  
}

ExpHGausVar <- function(sd, mean, K = 1){
  
  K * sd^2
  
}


x1 <- replicate(10000, HGaus(K = 14, sd = 5, mean = 1))
mean(x1)
var(x1)

ExpHGausVar(sd = 5, mean = 1)

ExpHGausVar(sd = 5, mean = 1, K = 14)

