# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(terra)
library(latex2exp)
library(nimble)
library(ASMbook)
library(lme4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make an empty polygon (in this class we're just going to make squares, but
#                        feel free to get wild on your own time)
# ng: dimensions of grid
# min: minimum dimension
# max: maximum dimension
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ng <- 50
min <- 0
max <- ng
box <- st_sfc(st_polygon(
  list(cbind(c(0,ng,ng,0,0),
             c(0,0,ng,ng,0))))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a grid
# nc: total number of cells in grid
# g: a grid
# c: grid centroids
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nc <- ng^2
g <- st_make_grid(box, n = c(ng, ng), cellsize = max/ng)
str(g)
p <- st_centroid(g)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize the grid we just made... (using base plot)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(0,0,0,0), oma = c(1,1,1,5))
plot(box)
plot(g, add = T)
# points(c(25,25,49) ~ c(1,25,49))
head(g)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use a geospatial model to generate something similar to a GRF
# *note: this is a complicated multivariate normal geospatial model
# don't worry too much about how this gets constructed (unless you
# find it fascinating and want to dive into it), just know that it 
# produces some random noise and cells that are close together are 
# more similar than random
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# D: CREATE A DISTANCE MATRIX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D <- st_distance(p, p)
dim(D)
par(mfrow = c(1,1), mar = c(5,5,2,2))
hist(D)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C: USE A FUNCTION TO CREATE A COVARIANCE MATRIX
# generate multivariate normal spatial random effects to fill grid (truth is z)
# spatial random effects
# https://www2.stat.duke.edu/~cr173/Sta444_Sp17/slides/Lec20.pdf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spatial.sigma <- 1
sigma2 <- 0.5
r <- 0.25
Sigma = spatial.sigma * exp(-r * D) + sigma2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use Cholesky decomposition to speed up this line of code
# epsilon <- MASS::mvrnorm(1, rep(0, nc), Sigma)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L <- chol(Sigma)
dim(Sigma)
x <- rnorm(nc, 0, 1)
epsilon <-  as.numeric(x %*% L)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign values to z
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- st_as_sf(g)
z$z <- epsilon
plot(z, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot grid-level sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s <- sample(nc, 100, replace = F)

plot(z, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
points(p[s], pch = 19)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot sub-grid-level sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s <- sample(nc, 10, replace = F)
s <- rep(s, each = 10)
plot(z, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
points(st_jitter(p[s], amount = 0.5), pch = 19)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here the box is our extent
# a single grid cell is our grain
# scope is ratio between extent and grain
#
# we're not really doing psuedo-replication here... there's no experiment
# but our sampling won't be treated appropriately in one of the models ('pseudo')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nSims <- 500

mu <- matrix(NA, nSims, 3)
se <- matrix(NA, nSims, 3)
diff <- matrix(NA, nSims, 3)
out <- matrix(NA, nSims, 3)

for (ii in 1:nSims){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # design A: 'normal'
  # sample (s) 100 cells
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n <- 100
  s <- sample(nc, n, replace = F)
  
  # measure cells with a bit of error
  y <- z$z[s] + rnorm(n, 0, 0.1)

  m <- lm(y ~ 1)
  # summary(m)
  # str(summary(m))
  mu[ii, 1] <- m$coefficients[1]
  se[ii, 1] <- summary(m)$coefficients[2] 
  diff[ii,1] <- mu[ii,1] - mean(z$z)
  out[ii,1] <- ifelse(mean(z$z) > mu[ii,1] + 1.96*se[ii,1] |
                        mean(z$z) < mu[ii,1] - 1.96*se[ii,1],
                      1, 0)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # design B: 'psuedo-replication'
  # there's no experiment here... but similar concept
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n <- 10
  s <- sample(nc, n, replace = F)
  s <- rep(s, each = 10)
  y <- z$z[s] + rnorm(length(s), 0, 0.1)
  
  m <- lm(y ~ 1)
  mu[ii, 2] <- m$coefficients[1]
  se[ii, 2] <- summary(m)$coefficients[2]   
  diff[ii,2] <- mu[ii,2] - mean(z$z)  
  out[ii,2] <- ifelse(mean(z$z) > mu[ii,2] + 1.96*se[ii,2] |
                        mean(z$z) < mu[ii,2] - 1.96*se[ii,2],
                      1, 0)  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # using random effects and the same data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  m <- lmer(y ~ 1 + (1|s))
  mu[ii, 3] <- summary(m)$coefficients[1] 
  se[ii, 3] <- summary(m)$coefficients[2] 
  diff[ii,3] <- mu[ii,3] - mean(z$z)  
  out[ii,3] <- ifelse(mean(z$z) > mu[ii,3] + 1.96*se[ii,3] |
                        mean(z$z) < mu[ii,3] - 1.96*se[ii,3],
                      1, 0) 
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How often is the true value outside the confidence intervals
# of the estimate?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colMeans(out)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What do estimates of standard error look like for the three models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5,5,2,2))
boxplot(se, names = c('n = 100','Pseduo','Random fx'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What do estimates of means look like for the three models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5,5,2,2))
boxplot(mu, names = c('n = 100','Pseduo','Random fx'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) in each instance, we're taking 100 samples... or are we?
# 2) why are the SEs greater for the random effects model?
# 3) why are the means the same for the 'pseudo' and 'random fx' model?
# 4) why do the 'normal' and 'random' effects model contain the true value
#    the same amount of the time, and why is the 'pseduo' approach so bad?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's do an actual experiment, we'll treat half of the sampled cells
# 50 or 5 depending on our sampling design...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nSims <- 500
b <- 1
beta <- matrix(NA, nSims, 3)
se <- matrix(NA, nSims, 3)
diff <- matrix(NA, nSims, 3)
out <- matrix(NA, nSims, 3)

for (ii in 1:nSims){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # design A: 'normal'
  # sample (s) 100 cells
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n <- 100
  s <- sample(nc, n, replace = F)
  t <- sample(s, n/2, replace = F)
  
  # measure cells with a bit of error
  # and add treatment effect
  y <- z$z[s] + rnorm(n, 0, 0.1)
  y[which(s %in% t)] <- y[which(s %in% t)] + b
  
  # create a covariate
  e <- rep(0, length(s))
  e[which(s %in% t)] <- 1
  
  
  m <- lm(y ~ e)
  # summary(m)
  # str(summary(m))
  beta[ii, 1] <- m$coefficients[2]
  se[ii, 1] <- summary(m)$coefficients[2,2] 
  diff[ii,1] <- beta[ii,1] - mean(z$z)
  out[ii,1] <- ifelse(b > beta[ii,1] + 1.96*se[ii,1] |
                        b < beta[ii,1] - 1.96*se[ii,1],
                      1, 0)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # design B: 'psuedo-replication'
  # there's no experiment here... but similar concept
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n <- 10
  s <- sample(nc, n, replace = F)
  t <- sample(s, n/2, replace = F)
  s <- rep(s, each = 10)
  y <- z$z[s] + rnorm(length(s), 0, 0.1)
  y[which(s %in% t)] <- y[which(s %in% t)] + b
  
  e <- rep(0, length(s))
  e[which(s %in% t)] <- 1
  
  m <- lm(y ~ e)
  beta[ii, 2] <- m$coefficients[2]
  se[ii, 2] <- summary(m)$coefficients[2,2]   
  diff[ii,2] <- beta[ii,2] - mean(z$z)  
  out[ii,2] <- ifelse(b > beta[ii,2] + 1.96*se[ii,2] |
                        b < beta[ii,2] - 1.96*se[ii,2],
                      1, 0)  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # using random effects and the same data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  m <- lmer(y ~ e + (1|s))
  beta[ii, 3] <- summary(m)$coefficients[2] 
  se[ii, 3] <- summary(m)$coefficients[2,2] 
  diff[ii,3] <- beta[ii,3] - mean(z$z)  
  out[ii,3] <- ifelse(b > beta[ii,3] + 1.96*se[ii,3] |
                        b < beta[ii,3] - 1.96*se[ii,3],
                      1, 0) 
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's plot the estiamtes of the treatment effect from every simulation...
# what are the key differences?
# why are they occurring?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,3))
plot(beta[,1], ylim = c(-2,5))
arrows(1:nSims, beta[,1]+1.96*se[,1], 1:nSims, beta[,1]-1.96*se[,1],length = 0)

plot(beta[,2], ylim = c(-2,5))
arrows(1:nSims, beta[,2]+1.96*se[,2], 1:nSims, beta[,2]-1.96*se[,2],length = 0)

plot(beta[,3], ylim = c(-2,5))
arrows(1:nSims, beta[,3]+1.96*se[,3], 1:nSims, beta[,3]-1.96*se[,3],length = 0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how often is the true effect (b = 1) outside the confidence intervals
# of the estimate of the effect with the three methods? Why?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colMeans(out)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# type I error: false positive, or erroneous rejection of true null hypothesis
#               here our null would be 'no effect of treatment'
#               go back up in the code, change b to 0, and assess the type I error
#               rate using colMeans(out)
#               
#               which method is most susceptible to type I error?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# type II error: false negative, or failure to reject a false null...
#                how might we assess this for each of the three models?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

