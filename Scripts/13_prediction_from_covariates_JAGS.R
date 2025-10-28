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
# contents
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) simulate a remotely-sensed landscape {Lines 20-100}
# (2) simulate abundance as a function of landscape (Lines 100-120)
# (3) sample abundance (Lines 100-120)
# (4) run models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) simulate a landscape {Lines 20-100}
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# D: CREATE A DISTANCE MATRIX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D <- st_distance(p, p)
dim(D)
par(mfrow = c(1,1), mar = c(5,5,2,2))
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
z$c <- epsilon
plot(z['c'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (2) simulate abundance as a function of landscape (Lines 100-130)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(1, 0.5)
E <- exp(alpha[1] + alpha[2] * epsilon + rnorm(nc, 0, 0.1))
y <- rpois(nc, E)
z$y <- y
plot(z["y"], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,20,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample the landscape (counts)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100
s <- sample(nc, 100, replace = F)

plot(z['c'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
points(p[s], pch = 19)


data <- data.frame(matrix(NA, nrow = n, ncol = 2))
names(data) <- c('y','x')
data$y <- z$y[s]
data$x <- z$c[s]

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(y ~ x, data = data, ylab = 'Abundance', xlab = 'Covariate', cex.lab = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# what we actually see from our abundance data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# z$y[-c(s)] <- NA
z$o <- z$y
z$o[-c(s)] <- NA
plot(z['o'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (4) run a model
# a) Likelihood-based model
# b) Bayesian model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a) ML estimates and predictions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(mL <- glm(y ~ x, data = data, family = 'poisson'))
pL <- predict(mL, newdata = (data.frame(x = z$c)), interval = 'confidence', se.fit = T)
z$pL <- exp(pL[[1]])

plot(z['pL'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,20,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)   


par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(z$y ~ z$pL, ylab = 'Abundance', xlab = 'predicted abundance')
abline(0,1)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# b) Bayesian estimates and predictions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("poisson.jags")
cat("
      model {

      alpha[1] ~ dnorm(2, 1)
      alpha[2] ~ dnorm(0, 1)

      for (i in 1:n){
        log(lambda[i]) = alpha[1] + alpha[2] * x[i]
        y[i] ~ dpois(lambda[i])
        
      }
      
      }
      ",fill = TRUE)
sink()

# data to pass to model
jags.data <- list(y = data$y, x = data$x, n = n)

# Initial values
inits <- function(){list(alpha = c(2,0.1))}  

# Parameters monitored
parameters <- c('alpha')

# MCMC settings
nc <- 4
nt <- 10
ni <- 10000
nb <- 5000


library(jagsUI)
Sys.time()
mB <- jags(jags.data, inits, parameters, "poisson.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()
summary(mB)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# output predictions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(z)
iters <- mB$mcmc.info$n.samples
px <- z$c
py <- matrix(NA, iters, ng^2)

for (j in 1:ncol(py)){
  py[,j] <- exp(mB$sims.list$alpha[,1] + mB$sims.list$alpha[,2] * px[j])
}
z$pB <- colMeans(py)
plot(z['pB'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,20,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)   


par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(z$pB ~ z$pL, ylab = 'Bayes predictions', xlab = 'ML predictions', cex.lab = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Think carefully about our priors. Make a plot of what our prior for 
#    the intercept looks like on the real scale (-Inf, Inf) and its exponent (0, Inf)
# 2) Try to provide informative enough priors that you actually change
#    the Bayesian predictions?
# 3) If you're feeling adventurous...
#    What is missing from our model that is present in the simulation?
#    Would it be possible to add that in?
#    Take a shot at it!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

