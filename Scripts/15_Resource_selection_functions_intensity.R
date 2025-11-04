# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(terra)
library(latex2exp)
library(nimble)
library(ASMbook)
library(lme4)
library(MuMIn)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# contents
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) simulate a remotely-sensed landscape {Lines 20-100}
# (2) simulate abundance as a function of landscape (Lines 100-120)
# (3) sample abundance (Lines 100-120)
# (4) run four models and calculate AICc, think about model weights
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
ndvi <-  as.numeric(x %*% L)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign values to z
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- st_as_sf(g)
z$ndvi <- ndvi
plot(z['ndvi'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a patchy cheatgrass invasion
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
b <- rnorm(nc, 0, 1)
risk <- as.numeric(b %*% L)
cheat <- rbinom(nc, 1, plogis(risk))
z$cheat <- cheat
plot(z['cheat'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,1,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (2) simulate abundance as a function of landscape (Lines 100-130)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(-1, 1, -1)
intensity <- exp(alpha[1] + alpha[2] * ndvi + alpha[3] * cheat + rnorm(nc, 0, 0.1))
y <- rpois(nc, intensity)
z$y <- y
plot(z["y"], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,30,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's plot our true outcomes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(jitter(y) ~ ndvi, ylab = 'Use', xlab = 'NDVI',
     col = cheat + 1, cex.lab = 2, pch = 19)

boxplot(y ~ cheat, col = c('grey', 'firebrick2'), names = c('Native', 'Cheatgrass'),
        xlab = '', ylab = 'Abundance (n)', cex.lab = 2)








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run a model in Nimble
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  for (i in 1:n){
    mu[i] <- exp(alpha[1] + alpha[2] * ndvi[i] + alpha[3] * cheat[i])
    y[i] ~ dpois(mu[i])
  }
  
  # priors
  alpha[1] ~ dnorm(0, 1)
  alpha[2] ~ dnorm(0, 0.1)
  alpha[3] ~ dnorm(0, 0.1)  
  
} )


str(dataList <- list(y = y, ndvi = ndvi, cheat = cheat, 
                     n = 2500)) 

inits <- function()
  list (alpha = c(1,0,0))

params <- c("alpha")

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m <- nimbleMCMC(code = model, constants = dataList,
                   inits = inits, monitors = params, nburnin = nb, niter = ni,
                   thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m, head) # print first 5 values of each chain (not shown)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) What are some key weaknesses of this approach?
#    a) where did we mark grouse? why does that matter?
#    b) how many individuals are marked? at what temporal scale are we collecting
#       data on marked individuals?
#    c) how could we address those issues?! :)
# 
# 2) In the code below, we extract posterior distributions and create two figures
#    How could we modify the code to create predictions for a subset of figures?
#    Specifically, imagine NDVI values of -2 and 2, and with and without
#    invasive cheatgrass... Solution (far) below
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USE OUR POSTERIOR DISTRIBUTIONS TO:
# predict intensity of use across grid cells
# predict intensity of use across a user-defined range of covariate values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- do.call(rbind.data.frame, m$samples)
iters <- nrow(tmp)
summary(ndvi)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a prediction for each grid cell
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
py <- matrix(NA, iters, ng^2)
for (j in 1:ncol(py)){
  py[,j] <- exp(tmp$'alpha[1]' + tmp$'alpha[2]' * ndvi[j] + tmp$'alpha[3]' * cheat[j])
}
z$p <- colMeans(py)
plot(z['p'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,20,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)   

# predictions vs. observations
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(z$y ~ z$p, ylab = 'Observations', xlab = 'Bayesian predictions', cex.lab = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derive expected cell values across a range of 
# covariate values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
resolution = 100
pndvi <- seq(-3, 3, length.out = resolution)
Ey <- array(NA, dim = c(iters, resolution,2))
qEy <- array(NA, dim = c(5, resolution, 2))
for (j in 1:resolution){
  Ey[,j,1] <- exp(tmp$'alpha[1]' + tmp$'alpha[2]' * pndvi[j] + tmp$'alpha[3]' * 0)
  Ey[,j,2] <- exp(tmp$'alpha[1]' + tmp$'alpha[2]' * pndvi[j] + tmp$'alpha[3]' * 1)  
  for (k in 1:2){
    qEy[,j,k] <- quantile(Ey[,j,k], c(0.025,0.25,0.5,0.75,0.975))
  }
}

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(qEy[3,,1] ~ pndvi, type = 'l', lwd = 2, ylim = c(0,10),
     ylab = "Intensity of use", xlab = 'NDVI', cex.lab = 1.5,
     las = 1)
lines(qEy[1,,1] ~ pndvi, lty = 2)
lines(qEy[5,,1] ~ pndvi, lty = 2)

lines(qEy[3,,2] ~ pndvi, lty = 1, col = 'red')
lines(qEy[1,,2] ~ pndvi, lty = 2, col = 'red')
lines(qEy[5,,2] ~ pndvi, lty = 2, col = 'red')

legend('topleft', legend = c('Native','Cheatgrass dominated'),
       lty = c(1,1), col = c('black','red'), bty = 'n', cex = 1.25)















# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SOLUTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
resolution = 3
pndvi <- seq(-3, 3, length.out = resolution)
Ey <- array(NA, dim = c(iters, resolution,2))
qEy <- array(NA, dim = c(5, resolution, 2))
for (j in 1:resolution){
  Ey[,j,1] <- exp(tmp$'alpha[1]' + tmp$'alpha[2]' * pndvi[j] + tmp$'alpha[3]' * 0)
  Ey[,j,2] <- exp(tmp$'alpha[1]' + tmp$'alpha[2]' * pndvi[j] + tmp$'alpha[3]' * 1)  
  for (k in 1:2){
    qEy[,j,k] <- quantile(Ey[,j,k], c(0.025,0.25,0.5,0.75,0.975))
  }
}

lx <- c(pndvi-0.1)
rx <- c(pndvi+0.1)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(qEy[3,,1] ~ lx, ylim = c(0,10), xlim = c(-4,4), xaxt = 'n',
     ylab = "Intensity of use", xlab = 'NDVI', cex.lab = 1.5,
     las = 1)
# add axis
axis(side = 1, at = pndvi, labels = pndvi)
# add CrIs for native
arrows(lx, qEy[1,,1], lx, qEy[5,,1], length = 0)
# add CrIs for cheatgrass-dominated
arrows(rx, qEy[1,,2], rx, qEy[5,,2], length = 0)
# add points
points(qEy[3,,1] ~ lx, pch = 21, bg = 'grey50', cex = 2)
points(qEy[3,,2] ~ rx, pch = 21, bg = 'red', cex = 2)

legend('topleft', legend = c('Native','Cheatgrass dominated'),
       pch = 21, pt.bg = c('grey50','red'), bty = 'n', cex = 1.25)







