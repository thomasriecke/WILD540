library(reshape2)
library(vioplot)
library(nimble)
library(ASMbook)

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a simple linear model, lm()
# with raw length data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(m1 <- lm(y ~ x))

plot(y ~ x, ylab = 'Mass (lbs)', xlab = 'Length (ft)', 
     xlim = c(-1,10), ylim = c(-100,200),
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
abline(m1$coefficients[1], m1$coefficients[2], lwd = 2, col = 'black')
points(0, m1$coefficients[1], pch = 21, bg = 'red', cex = 3) # the intercept

# a second point at a centered intercept
mean(x)
points(mean(x), m1$coefficients[1] + m1$coefficients[2] * mean(x), 
       pch = 21, bg = 'dodgerblue4', cex = 4)


# this code just shows the best fit if the intercept was fixed to 0...
# plot(y ~ x, ylab = 'Mass (lbs)', xlab = 'Length (ft)', 
#      xlim = c(-1,10), ylim = c(-100,200),
#      pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
# summary(m1b <- lm(y ~ 0 + x))
# abline(0, m1b$coefficients[1], lwd = 2, col = 'black')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll center x...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c <- x - mean(x)

plot(y ~ c, ylab = 'Mass (lbs)', xlab = 'Distance from average length (ft)', 
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)

summary(m2 <- lm(y ~ c))
# the intercept changed, but the slope is the same? why?

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll transform x (length in feet) to length in mm (m)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m <- c * 25.4 * 12

plot(y ~ m, ylab = 'Mass (lbs)', xlab = 'Distance from average length (mm)', 
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
summary(m3 <- lm(y ~ m))
# note that the intercept is the same, but the slope changed dramatically. Why?

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll z-standardize
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- (x - mean(x))/sd(x)
plot(y ~ z, ylab = 'Mass (lbs)', xlab = 'Z-standardized length', 
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
summary(m4 <- lm(y ~ z))
# The intercept is the same, the slope is slighlty different (it's now 
# on the scale of standard deviations!)











# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's write a Nimble model (or two...)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# so... now we have the same covariate on four different scales:
# x: alligator measurements in feet
# c: 'centered' alligator measurements in feet
# m: 'centered' alligator measurements in millimeters
# z: 'z-standardized' alligator measurements in feet
#
# We have a measurement of a response variable (y; mass in lbs)
# 
# and we have a sample size (n = 200)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean(y)


hist(rnorm(10000, 40, 10), breaks = 50, 
     xlab = 'Normal(40, sd = 10)',
     main = 'Prior for mass of an average alligator (lbs)')

hist(rgamma(10000, 1, 1), breaks = 50, 
     xlab = 'Gamma(1, 1)',
     main = 'Prior for error (lbs)')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the first step is writing the model itself, to do this we'll specify
# NimbleCode in our environment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model1 <- nimbleCode( {
  
  # Priors 
  intercept ~ dnorm(40, sd = 10)
  slope ~ dnorm(0, sd = 10)
  error ~ dgamma(1,1)
  # Likelihood
  for(i in 1:n){
    # mu: expected mass (lbs)
    mu[i] <- intercept + slope * z[i]
    # model
    y[i] ~ dnorm(mu[i], sd = error)
  }
} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = y, n = n, z = z)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list (intercept = 43, slope = 10, error = 2)

# Parameters monitored: same as before
params <- c("intercept", "slope", "error")

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m5 <- nimbleMCMC(code = model1, constants = dataList,
                       inits = inits, monitors = params, nburnin = nb, niter = ni,
                       thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m5, head) # print first 5 values of each chain (not shown)

par(mfrow = c(1, 3));
coda::traceplot(m5$samples) # not shown
m5$summary$all.chains # Posterior summaries from all 4 chains

m5sum <- nimble_summary(m5$samples, params) # Summary table
round(m5sum, 3)

# pretty close to frequentist estimates
summary(m4)
round(m5sum, 3)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's write the same model with greek letters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model2 <- nimbleCode( {
  
  # Priors 
  alpha ~ dnorm(40, sd = 10)
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
str(dataList <- list(y = y, n = n, z = z)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list (alpha = 43, beta = 10, sigma = 2)

# Parameters monitored: same as before
params <- c("alpha", "beta", "sigma")

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m6 <- nimbleMCMC(code = model2, constants = dataList,
                   inits = inits, monitors = params, nburnin = nb, niter = ni,
                   thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m5, head) # print first 5 values of each chain (not shown)

par(mfrow = c(1, 3));
coda::traceplot(m5$samples) # not shown
m5$summary$all.chains # Posterior summaries from all 4 chains

m6sum <- nimble_summary(m6$samples, params) # Summary table
round(m6sum, 3)

# pretty close to frequentist estimates
summary(m4)
round(m6sum, 3)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# things to try on your own time...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 1) replace the z-standardized values of length (z) with centered values (c), centered
#    values in millimeters (m), or real lengths (x). 
#    a) how do you expect the parameter estimates to changes?
#    b) should the priors change? how?
#
# 2) try some different priors... heck, try some really BAD priors :)
#    a) what happened to your posterior estimates
#    b) how bad do you need to make your priors to change inference
#
# 3) change the sample size (up or down or both)
#    a) how does that change your answer to question 2 above
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




