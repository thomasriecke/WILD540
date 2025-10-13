# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# An introduction to categorical covariates, 
# in which we explore fixed and random effects, play with 'Bergmann's rule' a bit,
# and think about how random effects behave
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(1234)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'serif')
library(jagsUI)
library(vioplot)
library(lme4)
library(nimble)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# randomly sample 'l' locations, 
# note to fix: pull the actual Red fox home range and generate realistic coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
l <- 30
# spatial coordinates (unused at at the moment, but may rework to demonstrate
# spatially correlated random effects)
s <- matrix(NA, l, 2) # number of rows = l, lat+long = 2 columns
s[,1] <- runif(l,-120,-80)
s[,2] <- runif(l, 35, 55)

# let's extract latitude from the spatial coordinates (x)
x <- as.numeric(scale(s[,2]))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# randomly sample 'n' foxes from each population and assign a group id
# and latitude to each fox...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lambda will represent the average number of samples
lambda <- 7
# we'll simulate the actual number of samples from a poisson distribution
n <- rpois(l, lambda)
n

# now let's generate the average size of an adult female red fox (mu.star)
# as well as the variance among populations (sigma[1])
mu.star <- 12
sigma <- NULL
sigma[1] <- 0.5

# now we'll simulate each populations mean as a function of mu.star, 
# sigma[1], and an effect of latitude on population size (beta)
beta <- 0
mu <- rnorm(l, mu.star + beta * x, sigma[1])

# let's plot these population level means as a function of z-standardized latitude
plot(mu ~ x, xlab = 'Z-standardized latitude', 
     ylab = 'Population mean mass (lbs)', cex.lab = 2)

# let's plot these population level means as a function of actual latitude
plot(mu ~ s[,2], xlab = "Latitude (N\u00B0)", 
     ylab = 'Population mean mass (lbs)', cex.lab = 2)


# simulate a population identifier for each fox
sum(n)
p <- rep(seq(1:l), times = n)
p

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's simulate each foxes mass
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigma[2] <- 0.5
y <- rnorm(length(p), mu[p], sigma[2])

# make a boxplot of fox mass as a function of population
boxplot(y ~ p, ylab = 'Female fox mass (lbs)', xlab = 'Population id', cex.lab = 2)
p
# regress fox mass against latitude
plot(y ~ s[,2][p], xlab = "Latitude (N\u00B0)",
     ylab = 'Female fox mass (lbs)', pch = 21, bg = 'red2', cex.lab = 2, cex = 1.5)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's run two ML models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(m1 <- lm(y ~ as.factor(p)))       # fixed effects
summary(m2 <- lmer(y ~ (1|as.factor(p)))) # random effects















# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  # model each populations' mean as a function of a grand.mean (mu.star) 
  # and variance among groups (sigma[1]^2)
  for (j in 1:(l+1)){
    mu[j] ~ dnorm(mu.star, tau[1])
  }
  
  # estimate each fox's mass as drawn from a random group intercept
  for (i in 1:I){
    y[i] ~ dnorm(mu[p[i]], tau[2])
  }
  
  # priors
  mu.star ~ dnorm(12, 0.1)
  sigma[1] ~ dgamma(1,1)
  sigma[2] ~ dgamma(1,1)
  tau[1] <- 1/(sigma[1] * sigma[1])
  tau[2] <- 1/(sigma[2] * sigma[2])
  
  
} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = y, l = l, 
                     I = sum(n), p = p)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list (sigma = c(2,1), mu.star = 12)

# Parameters monitored: same as before
params <- c("mu", "mu.star", "sigma",'tau')

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m3 <- nimbleMCMC(code = model, constants = dataList,
                  inits = inits, monitors = params, nburnin = nb, niter = ni,
                  thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m3, head) # print first 5 values of each chain (not shown)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's make a plot from our nimble output
# first we'll put all the chains together into one big chain
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- do.call(rbind.data.frame, m3$samples)

str(m3$summary$all.chains)
# summary statistics
m3$summary$all.chains


plot(m3$summary$all.chains[1:31,2], ylim = c(7,17), xlim = c(0,32),
     ylab = 'Mean mass (lbs)',
     xlab = 'Population', cex.lab = 2)
arrows(1:(l+1), m3$summary$all.chains[1:31,4], 
       1:(l+1), m3$summary$all.chains[1:31,5], length = 0, lty = 2, lwd = 1.5)
points(m3$summary$all.chains[1:31,2], pch = 21, bg = 'red2', cex = 2.5)


# note that the last chain is estimated without data 
# (i.e., it's a direct estimate of what a hypothetical 31st population would
# look like...)
points(y ~ p, pch = 1, cex = 0.5)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model2 <- nimbleCode( {
  
  # model each populations' mean as a function of a grand.mean (mu.star) 
  # and variance among groups (sigma[1]^2)
  for (j in 1:(l+1)){
    mu[j] ~ dgamma(1.2, 0.1)
  }
  
  # estimate each fox's mass as drawn from a random group intercept
  for (i in 1:I){
    y[i] ~ dnorm(mu[p[i]], tau)
  }
  
  # priors

  sigma ~ dgamma(1,1)
  tau <- 1/(sigma * sigma)

  
} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = y, l = l, 
                     I = sum(n), p = p)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list (sigma = c(2), mu.star = 12)

# Parameters monitored: same as before
params <- c("mu", "sigma",'tau')

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m4 <- nimbleMCMC(code = model2, constants = dataList,
                   inits = inits, monitors = params, nburnin = nb, niter = ni,
                   thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m4, head) # print first 5 values of each chain (not shown)




quantile(rgamma(100000,1.2,0.1), c(0.025,0.5,0.975))




# plot next to each other
par(mfrow = c(1,2))
plot(m3$summary$all.chains[1:31,2], ylim = c(7,17), xlim = c(0,32),
     ylab = 'Mean mass (lbs)',
     xlab = 'Population', cex.lab = 2)
arrows(1:(l+1), m3$summary$all.chains[1:31,4], 
       1:(l+1), m3$summary$all.chains[1:31,5], length = 0, lty = 2, lwd = 1.5)
points(m3$summary$all.chains[1:31,2], pch = 21, bg = 'red2', cex = 2.5)

plot(m4$summary$all.chains[1:31,2], ylim = c(7,17), xlim = c(0,32),
     ylab = 'Mean mass (lbs)',
     xlab = 'Population', cex.lab = 2)
arrows(1:(l+1), m4$summary$all.chains[1:31,4], 
       1:(l+1), m4$summary$all.chains[1:31,5], length = 0, lty = 2, lwd = 1.5)
points(m4$summary$all.chains[1:31,2], pch = 21, bg = 'red2', cex = 2.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# what's the difference?
# there's one obvious difference: the estimate for our hypothetical 31st population
# Do you see any other differences? See the plot immediately below for a hint
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m3$summary$all.chains[1:30,2] ~ m4$summary$all.chains[1:30,2], 
     ylab = 'Random effect estimates of pop. mass',
     xlab = 'Fixed effect estimates of pop. mass')
abline(0,1)















# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework 1:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# There's a slight signal here [in the last plot above], do you see it?
# The estimates are the same near the mean, but the residuals are slightly
# above the line at low weights, and below the line at high weights.
# 
# This is called 'shrinkage.' Put simply, the random effects are shrinking ever so
# slighlty towards the mean. Go back up to the beginning of the script where
# we define lambda, or the average number of foxes sampled per population. Change
# lambda to ~4 and run both models again. What happens? Why?
#
# Now change lambda to 35 (we'll measure more foxes) and run the 
# entire script again. There should be less shrinkage. What do you predict would
# happen if we measured ~1000 foxes in each population (don't do this it will take
# hours!)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 5:
#
# Now let's imagine we want to estimate variation among populations 
# around a trend-line that we estimate as a function of latitude?
#
# How would we do that?
# 
# First, we can't with fixed effects because the estimates of independent population
# variation and the effect of latitude would be confounded (don't take my word
# for it, give it a shot and see what happens)
#
# We can do this with random effects. We will estimate the means of each population
# as a random draw around the trend-line estimated using an effect of latitude,
# i.e., random residuals.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("fox_random_latitude.jags")
cat("
    model {

    # model each populations' mean as a function of a grand.mean (mu.star) 
    # and variance among groups (sigma[1]^2)
    # 
    # note that this is the only thing we need to change, we simply add an 
    # effect of x (which we call beta) to affect each population level mean
    # and supply a prior for beta below.
    for (j in 1:l){
      mu[j] ~ dnorm(mu.star + beta * x[j], tau[1])
    }

    # estimate each fox's mass as drawn from a random group intercept
    for (i in 1:I){
      y[i] ~ dnorm(mu[p[i]], tau[2])
    }

    mu.star ~ dnorm(12, 0.1)
    beta ~ dnorm(0, 0.1)
    sigma[1] ~ dgamma(1,1)
    sigma[2] ~ dgamma(1,1)
    tau[1] = 1/(sigma[1] * sigma[1])
    tau[2] = 1/(sigma[2] * sigma[2])

    }
    ",fill = TRUE)
sink()

#############################################################################################
# Bundle data
#############################################################################################

# we've now added x (latitude of each population) as data
jags.data <- list(y = y, l = l, I = length(p), p = p, x = x)

inits <- function(){list()}  

# Parameters monitored
parameters <- c('sigma','mu','mu.star','beta')

nc <- 4
nt <- 25
ni <- 50000
nb <- 25000


# Call JAGS from R 
# 7s for 50k iterations with an Intel i9-10900 10-core processor with some other computing going on as well...
library(jagsUI)
Sys.time()
m5 <- jags(jags.data, inits, parameters, "fox_random_latitude.jags", parallel = T, 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Sys.time()

print(m5)


plot(m5$q50$mu ~ x, ylim = c(9,15), xlim = c(-3,3),
     ylab = 'Mean mass (lbs)',
     xlab = 'Truth', cex.lab = 2)
arrows(x, m5$q2.5$mu, x, m5$q97.5$mu, length = 0, lty = 2, lwd = 1.5)
points(m5$q50$mu ~ x, pch = 21, bg = 'red2', cex = 2.5)







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now let's compare estimates from this model to our random effects only model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m5$q50$mu ~ m3$q50$mu,
     ylab = 'Estimates from m(random + latitude)',
     xlab = 'Estimates from m(random)', cex.lab = 2)
abline(0,1)

boxplot(m3$sims.list$sigma[,1], m5$sims.list$sigma[,1])




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework Part II:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Why are the estimates of the amount of variance between groups different
#    between the two models?
#
# 2) Rerun the entire script, but change the number of sampled populations
#    to 4? What happens to our estimates of the grand mean and among-group
#    variance? Why?
plot(m3$sims.list$sigma[,1] ~ m3$sims.list$mu.star, cex.lab = 1,
     xlab = expression(mu*'*'), ylab = expression(sigma[1]))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



