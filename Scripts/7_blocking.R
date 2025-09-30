# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(terra)
library(latex2exp)
library(nimble)
library(ASMbook)

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
# make a new box (the right hand side of the figure)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
block <- st_sfc(st_polygon(
  list(cbind(c(25,50,50,25,25),
             c(0,0,50,50,0))))
)
plot(block, add = T, col = 'red')
plot(g, add = T)
# blk is the cells inside the red area
# out are the cells inside the white area
blk <- st_contains(block, g)[[1]]
`%!in%` = Negate(`%in%`)
out <- which(c(1:nc) %!in% blk)

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
# test if positive definite [should be ok]
# all(eigen(Sigma)$values > 0)


epsilon <- MASS::mvrnorm(1, rep(0, nc), Sigma)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# randomly sample 50 points in each area
# sample 25 for treatments (e.g., exclosures)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n1 <- n2 <- 50
s1 <- sample(out, n1, replace = F)
s2 <- sample(blk, n2, replace = F)
nt1 <- nt2 <- 25
t1 <- sample(s1, nt1, replace = F)
t2 <- sample(s2, nt2, replace = F)

cell <- c(s1, s2)
area <- c(rep(0,n1), rep(1,n2))
trt <- rep(0, n1 + n2)
trt[which(cell %in% t1)] <- 1
trt[which(cell %in% t2)] <- 1
# create a data.frame
dat <- data.frame(cell, area, trt)
# So we have 50 points in each area
table(dat$area, dat$trt)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign values (v) to grid cells (z)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- st_as_sf(g)
trt <- rep(0, nc)
ar <- rep(0, nc)

trt[c(t1,t2)] <- 1; table(trt)
ar[blk] <- 1; table(ar)
beta <- c(1,2)

v <- epsilon + beta[1]*ar + beta[2] * trt
z$v <- v

plot(z, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
points(p[s1], cex = 0.5)
points(p[t1], pch = 21, bg = 'darkblue', cex = 0.5)
points(p[s2], cex = 0.5)
points(p[t2], pch = 21, bg = 'darkblue', cex = 0.5)

dat$y <- z$v[dat$cell]
dat$noise <- epsilon[dat$cell]



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new, we'll make a boxplot plot of our estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5,5,2,2), oma = c(0,0,0,0))
boxplot(dat$y ~ dat$trt + dat$area, ylab = 'Data',
        names = c('Control', 'Treatment', 'Control', 'Treatment'),
        xlab = 'Group', col = c('white','white','red','red'), cex.lab = 2,
        las = 1)


boxplot(dat$y ~ dat$trt + dat$area)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(dat)
m <- lm(y ~ area + trt, data = dat)
summary(m)


# obtain predictions from 'new data'
newdata <- data.frame(area = c(0,0,1,1), trt = c(0,1,0,1))
p <- predict(m, newdata = newdata, interval = 'confidence')
# what's going on here?
p <- data.frame(p)
plot(p$fit, xaxt = 'n', ylim = c(-2,6), xlim = c(0.5,4.5),
     xlab = 'Group', ylab = 'Estimate', cex.lab = 2)
axis(side = 1, at = 1:4, labels = c('control','Treatment','control','Treatment'))


# add actual data points
points(dat$y[dat$area == 0 & dat$trt == 0] ~ jitter(rep(1,n1-nt1)), cex = 0.5)
points(dat$y[dat$area == 0 & dat$trt == 1] ~ jitter(rep(2,n1-nt1)), cex = 0.5)
points(dat$y[dat$area == 1 & dat$trt == 0] ~ jitter(rep(3,n1-nt1)), cex = 0.5)
points(dat$y[dat$area == 1 & dat$trt == 1] ~ jitter(rep(4,n1-nt1)), cex = 0.5)

# add confidence intervals
arrows(1:4, p$lwr, 1:4, p$upr, length = 0, lwd = 4)
# add means
points(p$fit, pch = 21, bg = c('white','white','red','red'), cex = 2.5)







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  # Priors 
  alpha ~ dnorm(0, sd = 10)
  beta[1] ~ dnorm(0, sd = 10)
  beta[2] ~ dnorm(0, sd = 10)
  sigma ~ dgamma(1,1)
  # Likelihood
  for(i in 1:n){
    # mu: expected value
    mu[i] <- alpha + beta[1] * ar[i] + beta[2] * tr[i]
    # data are a function of expected value and error
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = dat$y, n = nrow(dat), 
                     tr = dat$trt, ar = dat$area)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list (alpha = 0, beta = c(0,0), sigma = 2)

# Parameters monitored: same as before
params <- c("alpha", "beta", "sigma")

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m2 <- nimbleMCMC(code = model, constants = dataList,
                   inits = inits, monitors = params, nburnin = nb, niter = ni,
                   thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m2, head) # print first 5 values of each chain (not shown)

par(mfrow = c(1, 3));
coda::traceplot(m2$samples) # not shown
m2$summary$all.chains # Posterior summaries from all 4 chains

m2sum <- nimble_summary(m2$samples, params) # Summary table
round(m2sum, 3)


