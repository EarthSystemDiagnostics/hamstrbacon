

## Effect of hierarchical levels on variance of final level


HGammaVar <- function(K, shape, mean, adjust.shape = TRUE){
  
  AdjustShape <- function(shape, K){
    shape * K #+ ((K - 1)/(shape +1)) 
  }
  
  scale <- mean / shape
  
  if (adjust.shape){
    shape <- AdjustShape(shape, K)
  } else {
    shape = shape
  }
  
  
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
mn <- 10
K <- 20

x1 <- replicate(200^2, HGammaVar(K = K, shape = shp, mean = mn, adjust.shape = TRUE))
x1.m <- matrix(x1, ncol = 100)
c.means <- colMeans(x1.m)
hist(c.means)
mean(c.means)
var2 <- function(x) {
  mean((x - mean(x))^2)
}

var3 <- function(x) {
  mean(x^2) - mean(x)^2
}

x1.var <- apply(x1.m, 2, var)
#mean(x1.var)
hist(x1.var, 25, xlim = c(0, max(x1.var)))

abline(v = ExpVar(shape = shp, mean = mn), col = "Red")
abline(v = mean(x1.var), col = "Blue")
abline(v = median(x1.var), col = "Green")


#####

#x1 <- replicate(100^2, HGammaVar(K = 10, shape = shp, mean = mn, adjust.shape = TRUE))
g.pars <- MASS::fitdistr(x1, "gamma")
g.pars
f.gamma <- tibble(
  x = seq(-1, max(x1), length.out = 1000),
  fttd = dgamma(x, shape = g.pars$estimate[["shape"]], rate = g.pars$estimate[["rate"]]),
  true = dgamma(x, shape = shp, rate = shp / mn)
)

h1 <- hist(x1, 100, plot = FALSE)
hist(x1, 100, freq = FALSE, ylim = c(0, max(c(h1$density, f.gamma$true, f.gamma$fttd))))
lines(fttd~x, data = f.gamma, col = "Red")
lines(true~x, data = f.gamma, col = "Blue")

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

