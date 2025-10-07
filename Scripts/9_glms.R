# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# First, we'll use exactly the same simulation code we built to analyze
# the alligator mass ~ legnth regression relationship, but now we'll use a 
# log-normal distribution in our model
#
# Second, we'll use a logit-link to include a covariate for the probability
# of band-recovery probability using our duck banding example given a declining 
# number of duck hunters...
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First, let's simulate some alligator mass and length data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(reshape2)
library(vioplot)
library(nimble)

set.seed(1234)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# base plot settings
# (1x1 grid, 'times new roman-ish' font, large margins on x- and y-axes)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), family = 'serif', mar = c(5.1,5.1,2.1,2.1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate data for linear model demonstration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample size (n)
n <- 200

# simulate length (x)
x <- rlnorm(n, 1.75, 0.15)
hist(x, breaks = 50, main = NULL, xlab = 'Alligator length (ft)',
     col = 'forestgreen', border = 'forestgreen', cex.lab = 2)

# true regression parameters on the log scale
beta0 <- 1.6
beta1 <- 0.375

# simulate mass (y) as a function of length (x)
y <- rlnorm(n, beta0 + beta1 * x, 0.1)
plot(y ~ x, ylab = 'Mass (lbs)', xlab = 'Length (ft)',
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)

# z-standardize mass as a 
z <- as.numeric(scale(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run a model in glm()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(tmp <- glm(log(y) ~ z))
summary(z)
res <- 100
newdata <- data.frame(z = seq(-2,4, length.out = res))
p <- predict(tmp, newdata, se.fit = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# back-transform covariate values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bt <- newdata$z * sd(x) + mean(x)
plot(exp(p[[1]]) ~ bt, xlab = 'Length (ft)',
     ylab = 'Mass (lbs)', las = 1, type = 'l', ylim = c(0,250))
lines(exp(p[[1]] + 1.96*p[[2]]) ~ bt, lty = 2)
lines(exp(p[[1]] - 1.96*p[[2]]) ~ bt, lty = 2)
points(y ~ x)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run a model in JAGS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  # Priors 
  alpha ~ dnorm(0, sd = 10)
  beta ~ dnorm(0, sd = 10)
  sigma ~ dgamma(1,1)
  # Likelihood
  for(i in 1:n){
    # mu: expected mass (lbs)
    mu[i] <- alpha + beta * z[i]
    # model
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = log(y), n = n, z = z)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list (alpha = 43, beta = 10, sigma = 2)

# Parameters monitored: same as before
params <- c("alpha", "beta", "sigma")

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m <- nimbleMCMC(code = model, constants = dataList,
                   inits = inits, monitors = params, nburnin = nb, niter = ni,
                   thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m, head) # print first 5 values of each chain (not shown)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions (PART I)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) What does our prior for beta0 represent here? How would this look
#    on the scale of our data? i.e.,
hist(beta0.prior <- exp(rnorm(100000, 0, 316.2278)))
#    is this a good prior? why or why not? how much of our prior support is less than 1?, i.e.,
length(which(beta0.prior < 1))/length(beta0.prior)
#    what might a better prior look like? What does beta0 represent in your own words?
#
# First, let's imagine we do have a pretty good idea what an alligator of average (5.7 ft)
# length weighs, and that's about 45 lbs...
# we could back-transform this value to the real (-Inf, Inf) scale that our model 
# will work on quite easily
log(45)
# so there's that number! 3.806662
# 
# now let's imagine we want some uncertainty around that, perhaps anywhere from 
# 25 to 100 pounds as 95% credible intervals on the prior?
hist(rlnorm(10000, log(45), 1)) # too broad
hist(rlnorm(10000, log(45), 0.01)) # way too informative :)
hist(rlnorm(10000, log(45), 0.1))  # maybe too informative?
hist(rlnorm(10000, log(45), 0.25)) # not so bad
# we could even use what we know about the median and variance of a log-normal
# to try to pick values that would correspond to our desired prior?
# https://en.wikipedia.org/wiki/Log-normal_distribution
#
#
#
#
# 2) What about our prior for beta1? What does that represent, are all of these values
#    reasonable? How might we modify this prior to be more reflective of our 
#    expectations and prior knowledge rather than an incredibly vague statement
#    such as normal(0, sigma^2 = 100000)
#
# 3) There is code below to make a plot of our expectations of mass given 
#    different values of length. It's very similar to the plot we used in
#    week 3 when reviewing linear models. How (and why) does it differ?
#    Review the use of the melt function...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- do.call(rbind.data.frame, m$samples)
# resolution of the plot (number of columns in grid)
res <- 100
niter <- nrow(tmp)
py <- matrix(NA, niter, res)

# assess range of covariate values, and assign axis marks.
summary(x)
plotx <- 3:9
plotz <- (plotx - mean(x))/sd(x)

# now we make a vector across the range of our observed covariate values
pz <- seq(min(z),max(z),length.out = res)

# fill our prediction matrix for weight with uncertainty...
for (j in 1:res){
  py[,j] <- exp(m1$sims.list$beta0 + m1$sims.list$beta1 * pz[j])
}


# create a matrix to store credible (95% and 80%) intervals
qy <- matrix(NA, res, 5)
for (j in 1:res){
  qy[j,] <- quantile(py[,j], c(0.025,0.1,0.5,0.9,0.975))
}

# melt py
predy <- melt(py)
names(predy) <- c('iteration','pz','E')

# make a density plot of the expected relationship
# note you can adjust the colramp to inform how the densities are plotted
smoothScatter(predy$E ~ pz[predy$pz],                                     # y ~ and x coordinates
              colramp = colorRampPalette(c('white','forestgreen')), # color ramp
              nrpoints = 0,                                         # no little dark points
              xaxt = 'n',                                           # no x-axis
              ylab = 'Expected alligator mass (lbs)',
              xlab = 'Alligator length (ft)', cex.lab = 2)    
axis(side = 1, at = plotz, labels = plotx) # this is where our axis-labels come back into play

# add median
lines(qy[,3] ~ pz, lwd = 3, col = 'white', lty = 1)
# add 80% credible intervals
lines(qy[,2] ~ pz, lwd = 3, col = 'white', lty = 2)
lines(qy[,4] ~ pz, lwd = 3, col = 'white', lty = 2)
# add 95% credible intervals
lines(qy[,1] ~ pz, lwd = 3, col = 'white', lty = 3)
lines(qy[,5] ~ pz, lwd = 3, col = 'white', lty = 3)

# add actual data
points(y ~ z)












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Second, let's simulate some ducks being marked, released, and eventually
# shot, recovered, and reported by hunters...
# 
# We'll have a single covariate, the number of hunters. We might reasonably
# assume that there is a causal linkage between the number of hunters and 
# the proportion of ducks that are shot (note: if we don't at least test this 
# hypothesis we should probably just all throw in the towel on this whole 
# 'applied science' thing!) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# number of time periods, or years
nT <- 20
# number of releases each year [mark release 1000 ducks every year]
nR <- rep(1000, nT)
# number of hunters
nH <- rpois(nT, exp(7 + -0.125 * 1:nT))
plot(nH)

# z-standardize the number of hunters
z <- as.numeric(scale(nH))

# build logit-link model to estimate
beta0 <- -2
beta1 <- 0.25
f <- plogis(beta0 + beta1 * z)
plot(f ~ nH)

# simulate the number of actual recoveries and plot 
# number of recoveries divided by number of releases...
y <- rbinom(nT, nR, f)
plot(c(y/nR), ylab = 'ML band-recovery probability', xlab = 'Year', cex.lab = 2)




model <- nimbleCode( {
  
  # Priors 
  alpha ~ dnorm(0, sd = 100)
  beta ~ dnorm(0, sd = 100)
  # Likelihood
  for(i in 1:n){
    # mu: expected mass (lbs)
    logit(f[i]) <- alpha + beta * z[i]
    # model
    y[i] ~ dbin(f[i], N[i])
  }
} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = y, n = nT, z = z, N = nR)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list ()

# Parameters monitored: same as before
params <- c("alpha", "beta")

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m <- nimbleMCMC(code = model, constants = dataList,
                  inits = inits, monitors = params, nburnin = nb, niter = ni,
                  thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m, head) # print first 5 values of each chain (not shown)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework PART II
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) What does our prior for beta0 represent here? How would this look
#    on the scale of our data? i.e.,
hist(beta0.prior <- plogis(rnorm(100000, 0, 316.2278)), main = NULL, breaks = 100)
#    is this a good prior? why or why not? 
#    what might a better prior look like? What does beta0 represent in your own words?
#
# 2) What about our prior for beta1? What does that represent, are all of these values
#    reasonable? How might we modify this prior to be more reflective of our 
#    expectations and prior knowledge rather than an incredibly vague statement
#    such as normal(0, sigma^2 = 100000)
#
# 3) attempt to write code below to make a density plot 
#    along with credible intervals of expected values of band-recovery
#    probability (or the expected number of recoveries... or both?!) as a function
#    of the number of hunters
#
# 4) If you're really feeling wild, try to put a random effect on band-recovery
#    probability around the effect of the long-term decline in hunters?!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




