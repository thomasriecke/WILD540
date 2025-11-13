# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script simulates Resource Selection Function data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(nimble)
library(latex2exp)
library(sf)
set.seed(123)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# contents
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) simulate a remotely-sensed landscape {Lines 20-100}
# (2) simulate abundance as a function of landscape (Lines 100-120)
# (3) sample abundance (Lines 100-120)
# (4) run four models and calculate AICc, think about model weights
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
epsilon <- rnorm(nc, 0, 1)
ndvi <- as.numeric(epsilon %*% L)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign covariate values (X) to z
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- st_as_sf(g)
z$ndvi <- ndvi
colramp <- colorRampPalette(c('burlywood2','forestgreen'))
plot(z['ndvi'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate 2000 random points and extract covariate values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
random <- 2000
r <- st_sample(g, random)
points(r)
xyr <- st_intersects(r, g)
xyr <- unlist(xyr)

plot(z['ndvi'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)      
points(r, pch = 1, cex = 0.5)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate 500 used points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- 0; beta = 1
E <- plogis(alpha + beta * ndvi + rnorm(nc, 0, 0.1))
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(E ~ ndvi)

used <- 500
xyu <- NULL
u <- NULL
for (i in 1:500){
  # pick used cells
  xyu[i] <- which(rmultinom(1, 1, E/sum(E)) == 1)
  # randomly located used points within used cells
  u[i] <- st_sample(g[xyu[i]], 1)
}

u <- st_sfc(u)

plot(z['ndvi'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)      
points(u, pch = 19)
points(r, pch = 1, cex = 0.5)

ux <- ndvi[xyu]
rx <- ndvi[xyr]


par(mfrow = c(1,1), mar = c(5,5,2,2))
vioplot::vioplot(ux, rx, ndvi, drawRect = F,
                 names = c('Used','Random','Landscape'),
                 ylab = '')
mtext(TeX('Normalized Difference Vegetation Index (NDVI)'), 
      side = 2, line = 2.5, cex = 1.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create data....
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y <- c(rep(1, used), rep(0, random))
x <- c(ux, rx)
n = used + random

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# maximum likelihood model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(glm(y ~ x, family = 'binomial'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nimble model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  alpha ~ dnorm(0, sd = 10)
  beta ~ dnorm(0, sd = 10)
  for (i in 1:n){
    logit(theta[i]) <- alpha + beta * x[i]
    y[i] ~ dbern(theta[i])
  }

} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = y, x = x,  n= n)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()  list ()

# Parameters monitored: same as before
params <- c("alpha","beta")

# MCMC settings
ni <- 10000 ; nb <- 5000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m <- nimbleMCMC(code = model, constants = dataList,
                   inits = inits, monitors = params, nburnin = nb, niter = ni,
                   thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m, head) # print first 5 values of each chain (not shown)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samps <- do.call(rbind.data.frame, m$samples)
resolution <- 100
px <- seq(-3,3,length.out = resolution)
py <- matrix(NA, nrow(samps), resolution)
qy <- matrix(NA, 5, resolution)

for (j in 1:resolution){
  py[,j] <- plogis(samps$alpha + samps$beta * px[j])
  qy[,j] <- quantile(py[,j], c(0.025,0.25,0.5,0.75,0.975))
}

shape.type = c(1,19)

plot(jitter(y, 0.15) ~ x, cex = 0.5, cex.lab = 1.25,
     ylab = 'Resource selection function (RSF)', 
     xlab = 'Normalized Difference Vegetation Index (NDVI)',
     pch = shape.type[y+1])
lines(qy[3,] ~ px, lwd = 3)
lines(qy[1,] ~ px, lwd = 2, lty = 2)
lines(qy[5,] ~ px, lwd = 2, lty = 2)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here is code to write our own MCMC sampler
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0a) function to calculate log-likelihood
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LL <- function(param, y, x) {
  alpha <- param[1] # Intercept
  beta <- param[2] # Slope
  theta <- plogis(alpha + beta * x) # expected value of log-linear model
  LLi <- y * log(theta) + (1 - y) * log(1 - theta)
  LL <- sum(LLi) # LL for all observations in the data vector y
  return(LL)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate log-likelihood for two lines (MLE and generic 'bad')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(glm(y ~ x, family = 'binomial'))
LL(c(-1.5,0.35), y, x)
LL(c(0.5,-0.2), y, x)
LL(c(-1,-1), y, x)
px <- seq(-3,3, length.out = 100)


plot(jitter(y, 0.15) ~ x, cex = 0.5, cex.lab = 1.25,
     ylab = 'Resource selection function (RSF)', 
     xlab = 'Normalized Difference Vegetation Index (NDVI)',
     pch = shape.type[y+1])
lines(plogis(-1.5 + 0.35*px) ~ px, lty = 1, lwd = 3, col = 'black')
lines(plogis(0.5 + -0.2*px) ~ px, lty = 1, lwd = 3, col = 'red')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0b) function to calculate un-normalized log posterior density
# (we'll ignore the denominator, or probability of the data, as it doesn't vary)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.posterior <- function(alpha, beta, y, x){
  loglike <- LL(c(alpha, beta), y, x)
  logprior.alpha <- dnorm(alpha, mean = 0, sd = 10, log = TRUE)
  logprior.beta <- dnorm(beta, mean = 0, sd = 10, log = TRUE)
  return(loglike + logprior.alpha + logprior.beta)
}

# initial log posterior
log.posterior(-1,-1,y,x)
par(mfrow = c(1,1))
plot(jitter(y, 0.15) ~ x, cex = 0.5, cex.lab = 1.25,
     ylab = 'Resource selection function (RSF)', 
     xlab = 'Normalized Difference Vegetation Index (NDVI)',
     pch = shape.type[y+1])
lines(plogis(-1 + -1*px) ~ px, lty = 1, lwd = 3, col = 'black')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Initialize chains for both parameters: these are the initial values
# for alpha and beta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- -1
beta <- -1 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Choose number of iterations for which to run the algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
niter <- 10000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Create an R object to hold the posterior draws produced
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out.mcmc <- matrix(NA, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) create an R object to hold acceptance indicator for both params
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
acc <- rep(0, niter)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4b) create R objects to store every candidate proposal and log prior
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all <- matrix(NA, niter, 2)
lpr <- matrix(NA, niter, 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Evaluate current value of the log(posterior density) (for the inits)
# and plot current expected value and data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
logpost.curr <- log.posterior(alpha, beta, y, x)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) Choose values for the tuning parameters of the algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sd.tune <- c(0.1, 0.1) # first is for alpha, second for beta
logpost.cand <- NULL



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# demonstrate a first proposal
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a1 <- -1
b1 <- -1
ll1 <- LL(c(a1,b1), y, x)
a2 <- rnorm(1, a1, sd.tune[1])
b2 <- rnorm(1, b1, sd.tune[1])
ll2 <- LL(c(a2,b2), y, x)
lpa <- dnorm(alpha, mean = 0, sd = 10, log = TRUE)
lpb <- dnorm(alpha, mean = 0, sd = 10, log = TRUE)
ll2 + lpa + lpb

par(mfrow = c(1,1), mar = c(5,5,2,2), oma = c(0,0,0,0))
plot(jitter(y, 0.15) ~ x, cex = 0.5, cex.lab = 1.5,
     ylab = 'Resource selection function (RSF)', 
     xlab = 'Normalized Difference Vegetation Index (NDVI)',
     pch = shape.type[y+1])
lines(plogis(a1 + b1*px) ~ px, lty = 1, lwd = 3, col = 'black')
lines(plogis(a2 + b2*px) ~ px, lty = 1, lwd = 3, col = 'grey40')


exp(ll2 - ll1)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run MCMC algorithm for 10k iterations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(t in 1:niter){
  if(t %% 1000 == 0) # counter
    cat("iter", t, "\n")
  # First, update log-linear intercept (alpha)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # note that here we'll use Random walk block sampling to propose
  # a new value for alpha and beta at the same time
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Propose candidate value of alpha and beta at the same time!
  alpha.cand <- rnorm(1, alpha, sd.tune[1]) # note tuning 'parameter'
  beta.cand <- rnorm(1, beta, sd.tune[2]) # note tuning 'parameter'
  all[t,] <- c(alpha.cand, beta.cand)
  # Evaluate log(posterior) for proposed new alpha and
  # for current beta for the data y and covs x
  logpost.cand[t] <- log.posterior(alpha.cand, beta.cand, y, x)
  # Compute Metropolis acceptance ratio r
  r <- exp(logpost.cand[t] - logpost.curr)
  # Keep candidate alpha if it meets criterion (u < r)
  if(runif(1) < r){
    alpha <- alpha.cand
    beta <- beta.cand
    logpost.curr <- logpost.cand[t]
    acc[t] <- 1 # Indicator for whether candidate alpha accepted
  }
  
   out.mcmc[t,] <- c(alpha, beta) # Save samples for iteration i
   lpr[t,1] <- dnorm(alpha, mean = 0, sd = 10, log = TRUE)
   lpr[t,2] <- dnorm(beta, mean = 0, sd = 10, log = TRUE)
  # NOTE: if proposed new values not accepted, then in this step
  # we will copy their 'old' values and insert them into object 'out'
}



# Compute overall acceptance ratio separately for alpha and beta
acc.ratio <- mean(acc) # Should be in range 25%–40%
# Check the traceplots for convergence of chains
# Plot all MC draws right from the start (Fig. 2.15)
# Repeat traceplots only for post-burnin draws (Fig. 2.16)
par(mfrow = c(2, 1))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot progress at different samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samp <- 20
par(mfrow = c(1, 3))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(out.mcmc[samp,1] ~ samp, col = 'red', cex = 3)
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(out.mcmc[samp,2] ~ samp, col = 'red', cex = 3)
plot(jitter(y, 0.15) ~ x, cex = 0.5, cex.lab = 1.5,
     ylab = 'Resource selection function (RSF)', 
     xlab = 'Normalized Difference Vegetation Index (NDVI)',
     pch = shape.type[y+1])
lines(plogis(out.mcmc[samp,1] + out.mcmc[samp,2]*px) ~ px, lty = 1, lwd = 3, col = 'red')
log.posterior(out.mcmc[samp,1],out.mcmc[samp,2], y, x)


samp <- 200
par(mfrow = c(1, 3))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(out.mcmc[samp,1] ~ samp, col = 'red', cex = 3)
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(out.mcmc[samp,2] ~ samp, col = 'red', cex = 3)
plot(jitter(y, 0.15) ~ x, cex = 0.5, cex.lab = 1.5,
     ylab = 'Resource selection function (RSF)', 
     xlab = 'Normalized Difference Vegetation Index (NDVI)',
     pch = shape.type[y+1])
lines(plogis(out.mcmc[samp,1] + out.mcmc[samp,2]*px) ~ px, lty = 1, lwd = 3, col = 'red')
log.posterior(out.mcmc[samp,1], out.mcmc[samp,2], y, x)

samp <- 2000
par(mfrow = c(1, 3))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(out.mcmc[samp,1] ~ samp, col = 'red', cex = 3)
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(out.mcmc[samp,2] ~ samp, col = 'red', cex = 3)
plot(jitter(y, 0.15) ~ x, cex = 0.5, cex.lab = 1.5,
     ylab = 'Resource selection function (RSF)', 
     xlab = 'Normalized Difference Vegetation Index (NDVI)',
     pch = shape.type[y+1])
lines(plogis(out.mcmc[samp,1] + out.mcmc[samp,2]*px) ~ px, lty = 1, lwd = 3, col = 'red')
log.posterior(out.mcmc[samp,1], out.mcmc[samp,2], y, x)







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions about posteriors...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# burn in
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1, 2))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

burn <- 1000
out.burn <- out.mcmc[-c(1:burn),]
niter <- niter - burn

plot(1:niter, out.burn[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
plot(1:niter, out.burn[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how do we interpret posterior distributions?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qa <- quantile(out.burn[,1], c(0.025,0.25,0.5,0.75,0.975))
qb <- quantile(out.burn[,2], c(0.025,0.25,0.5,0.75,0.975))
length(out.burn[,2] > 0)/niter

par(mfrow = c(1, 2))
vioplot::vioplot(out.burn[,1], names = '', drawRect = F)
arrows(1, qa[1], 1, qa[5], length = 0, lwd = 2)
arrows(1, qa[2], 1, qa[4], length = 0, lwd = 6)
points(qa[3] ~ 1, cex = 4, pch = 19)
vioplot::vioplot(out.burn[,2], names = '', drawRect = F)
arrows(1, qb[1], 1, qb[5], length = 0, lwd = 2)
arrows(1, qb[2], 1, qb[4], length = 0, lwd = 6)
points(qb[3] ~ 1, cex = 4, pch = 19)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log-prior
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1, 3))
plot(lpr[,1], ylab = TeX("l($\\alpha$)"), xlab = 'Iteration', cex.lab = 1.5, type = 'l')
plot(lpr[,2], ylab = TeX("l($\\beta$)"), xlab = 'Iteration', cex.lab = 1.5, type = 'l')
w <- seq(-5, 5, length.out = 100)
plot(dnorm(w, mean = 0, sd = 10, log = T) ~ w, 
     ylab = 'Log probability density: normal(0, sd = 10)',
     xlab = TeX("$\\theta$"), cex.lab = 1.5,
     type = 'l')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lp of different priors?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1, 3))
lpr.new <- lpr
lpr.new[,1] <- dnorm(out.mcmc[,1], mean = 0, sd = 1, log = T)
lpr.new[,2] <- dnorm(out.mcmc[,2], mean = 0, sd = 1, log = T)
plot(lpr.new[,1], ylab = TeX("l($\\alpha$)"), xlab = 'Iteration', cex.lab = 1.5, type = 'l')
plot(lpr.new[,2], ylab = TeX("l($\\beta$)"), xlab = 'Iteration', cex.lab = 1.5, type = 'l')
w <- seq(-5, 5, length.out = 100)
plot(dnorm(w, mean = 0, sd = 1, log = T) ~ w, 
     ylab = 'Log probability density: normal(0, sd = 10)',
     xlab = TeX("$\\theta$"), cex.lab = 1.5,
     type = 'l')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# look at log-likelihood and log prob. density of priors
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
loglik <- NULL
for (i in 1:(niter+burn)){
  loglik[i] <- LL(c(out.mcmc[i,1],out.mcmc[i,2]), y, x)
}

par(mfrow = c(1, 3))
plot(loglik, ylab = TeX("l($\\alpha,\\beta$; y)"), 
     xlab = 'Iteration', cex.lab = 1.5, type = 'l')
plot(lpr[,1], ylab = TeX("l($\\alpha$)"), xlab = 'Iteration', cex.lab = 1.5, type = 'l')
plot(lpr[,2], ylab = TeX("l($\\beta$)"), xlab = 'Iteration', cex.lab = 1.5, type = 'l')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate likelihood surface
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100
alphaX <- seq(-2.5, -0.5, length.out = res)
betaX <- seq(-1.25, 1, length.out = res)

logL <- matrix(NA, res, res)

for (i in 1:res){
  for (j in 1:res){
    logL[i,j] <- LL(c(alphaX[i], betaX[j]), y, x)
  }
}
par(mfrow = c(1,1), mar = c(5,5,2,2))
fields::image.plot(logL, xaxt = 'n', yaxt = 'n')
axis(side = 1, labels = seq(min(alphaX), max(alphaX), length.out = 5), 
     at = seq(0,1,length.out = 5))
axis(side = 2, labels = seq(min(betaX), max(betaX), length.out = 5), 
     at = seq(0,1,length.out = 5))
mtext(TeX("$\\beta$"), side = 2, line = 3, cex = 2)
mtext(TeX("$\\alpha$"), side = 1, line = 3, cex = 2)














# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run MCMC algorithm for 10k iterations with different tuning parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alpha <- -1
beta <- -1 
niter <- 10000
out.mcmc <- matrix(NA, niter, 2, dimnames = list(NULL, c("alpha", "beta")))
acc <- rep(0, niter)
all <- matrix(NA, niter, 2)
lpr <- matrix(NA, niter, 2)
logpost.curr <- log.posterior(alpha, beta, y, x)
logpost.cand <- NULL
sd.tune <- c(0.1, 0.1) # first is for alpha, second for beta

for(t in 1:niter){
  if(t %% 1000 == 0) # counter
    cat("iter", t, "\n")
  # First, update log-linear intercept (alpha)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # note that here we'll use Random walk block sampling to propose
  # a new value for alpha and beta at the same time
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Propose candidate value of alpha and beta at the same time!
  alpha.cand <- rnorm(1, alpha, sd.tune[1]) # note tuning 'parameter'
  beta.cand <- rnorm(1, beta, sd.tune[2]) # note tuning 'parameter'
  all[t,] <- c(alpha.cand, beta.cand)
  # Evaluate log(posterior) for proposed new alpha and
  # for current beta for the data y and covs x
  logpost.cand[t] <- log.posterior(alpha.cand, beta.cand, y, x)
  # Compute Metropolis acceptance ratio r
  r <- exp(logpost.cand[t] - logpost.curr)
  # Keep candidate alpha if it meets criterion (u < r)
  if(runif(1) < r){
    alpha <- alpha.cand
    beta <- beta.cand
    logpost.curr <- logpost.cand[t]
    acc[t] <- 1 # Indicator for whether candidate alpha accepted
  }
  
  out.mcmc[t,] <- c(alpha, beta) # Save samples for iteration i
  lpr[t,1] <- dnorm(alpha, mean = 0, sd = 10, log = TRUE)
  lpr[t,2] <- dnorm(beta, mean = 0, sd = 10, log = TRUE)
  # NOTE: if proposed new values not accepted, then in this step
  # we will copy their 'old' values and insert them into object 'out'
}



# Compute overall acceptance ratio separately for alpha and beta
acc.ratio <- mean(acc) # Should be in range 25%–40%
# Check the traceplots for convergence of chains
# Plot all MC draws right from the start (Fig. 2.15)
# Repeat traceplots only for post-burnin draws (Fig. 2.16)
par(mfrow = c(2, 1), mar = c(5,5,5,2))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

mean(acc)


plot(logpost.cand[acc == 1], ylim = c(-1250,-1235), type = 'l',
     ylab = 'Log posterior')
points(logpost.cand[acc == 0], type = 'l', add = T, col = adjustcolor('red', alpha = 0.5))




