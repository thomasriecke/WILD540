# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load a package associated with Marc Kery and Ken Kellner's book
# we'll use the beeeater data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ASMbook)
library(latex2exp)
set.seed(123)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Counts of bee-eaters in Switzerland and actual (a) and scaled (x)
# year covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y <- c(0, 2, 7, 5, 1, 2, 4, 8, 10, 11, 16, 11, 10, 13, 
       19, 31, 20, 26, 19, 21, 34, 35, 43, 53, 66, 61, 72, 
       120, 102, 159, 199)
year <- 1990:2020
x <- (year - 1989) - 16


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a figure of counts across time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5,5,2,2))
bee.cols <- colorRampPalette(c("turquoise",'yellow','deeppink4'))
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to calculate log-likelihood
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LL <- function(param, y, x) {
  alpha <- param[1] # Intercept
  beta <- param[2] # Slope
  lambda <- exp(alpha + beta * x) # expected value of log-linear model
  LLi <- dpois(y, lambda, log=TRUE) # LL contribution of each datum i
  LL <- sum(LLi) # LL for all observations in the data vector y
  return(LL)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate log-likelihood of two different parameter pairs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
good.params <- c(2.886821, 0.147230)
bad.params <- c(3, -0.1)
lines(exp(good.params[1] + good.params[2] * x) ~ year, lwd = 3)
lines(exp(bad.params[1] + bad.params[2] * x) ~ year, lwd = 3, col = 'red')

# calculate log-likelihoods for good parameters (ML estimates) and bad params
LL(good.params, y, x)
LL(bad.params, y, x)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to calculate un-normalized log posterior density
# (we'll ignore the denominator, or probability of the data, as it doesn't vary)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.posterior <- function(alpha, beta, y, x){
  loglike <- LL(c(alpha, beta), y, x)
  logprior.alpha <- dnorm(alpha, mean = 0, sd = 100, log = TRUE)
  logprior.beta <- dnorm(beta, mean = 0, sd = 100, log = TRUE)
  return(loglike + logprior.alpha + logprior.beta)
}


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
acc <- matrix(0, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Evaluate current value of the log(posterior density) (for the inits)
# and plot current expected value and data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
logpost.curr <- log.posterior(alpha, beta, y, x)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha + beta * x) ~ year, lwd = 3)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) Choose values for the tuning parameters of the algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sd.tune <- c(0.1, 0.01) # first is for alpha, second for beta


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a 3x1 plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(3,1), mar = c(5,5,2,2))
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha + beta * x) ~ year, lwd = 3)
plot(alpha, xlim = c(0,niter), ylab = TeX("$\\alpha$"), cex.lab = 2, ylim = c(-3.25,3.25))
abline(h = 2.886821, col = 'grey70', lwd = 2)
plot(beta, xlim = c(0,niter), ylab = TeX("$\\alpha$"), cex.lab = 2, ylim = c(-1,0.2))
abline(h = 0.147230, col = 'grey70', lwd = 2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7 & 8) perform a single additional iteration for alpha
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample alpha a single time and calculate Metropolic acceptance ratio (r)
alpha.cand <- rnorm(1, alpha, sd.tune[1])
logpost.cand <- log.posterior(alpha.cand, beta, y, x)
r <- exp(logpost.cand - logpost.curr)
if(runif(1) < r){
  alpha <- alpha.cand
  logpost.curr <- logpost.cand
}
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha + beta * x) ~ year, lwd = 3)
lines(exp(alpha.cand + beta * x) ~ year, lwd = 3, col = 'grey50')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9) sample beta a single time and calculate Metropolic acceptance ratio (r)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta.cand <- rnorm(1, beta, sd.tune[2])
logpost.cand <- log.posterior(alpha, beta.cand, y, x)

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha + beta * x) ~ year, lwd = 3)
lines(exp(alpha + beta.cand * x) ~ year, lwd = 3, col = 'grey50')

r <- exp(logpost.cand - logpost.curr)
# Keep candidate if it meets criterion (u < r)
if(runif(1) < r){
  beta <- beta.cand
  logpost.curr <- logpost.cand
}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# overwrite initial values and Run MCMC algorithm for 10k iterations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- -1
beta <- -1 

for(t in 1:niter){
  if(t %% 1000 == 0) # counter
    cat("iter", t, "\n")
  # First, update log-linear intercept (alpha)
  # ––––––––––––––––––––––––––––––––––––––––––
  # Propose candidate value of alpha
  alpha.cand <- rnorm(1, alpha, sd.tune[1]) # note tuning 'parameter'
  # Evaluate log(posterior) for proposed new alpha and
  # for current beta for the data y and covs x
  logpost.cand <- log.posterior(alpha.cand, beta, y, x)
  # Compute Metropolis acceptance ratio r
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate alpha if it meets criterion (u < r)
  if(runif(1) < r){
    alpha <- alpha.cand
    logpost.curr <- logpost.cand
    acc[t,1] <- 1 # Indicator for whether candidate alpha accepted
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Second, update log-linear slope (beta)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # propose a candidate beta-value
  beta.cand <- rnorm(1, beta, sd.tune[2])

  # Evaluate the log(posterior) for proposed new beta and
  # for the current alpha for the data y with covs x
  logpost.cand <- log.posterior(alpha, beta.cand, y, x)
  # Compute Metropolis acceptance ratio
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta <- beta.cand
    logpost.curr <- logpost.cand
    acc[t,2] <- 1 # Indicator for whether candidate beta accepted
  }
  out.mcmc[t,] <- c(alpha, beta) # Save samples for iteration i
  # NOTE: if proposed new values not accepted, then in this step
  # we will copy their 'old' values and insert them into object 'out'
}



# Compute overall acceptance ratio separately for alpha and beta
(acc.ratio <- apply(acc, 2, mean)) # Should be in range 25%–40%
# Check the traceplots for convergence of chains
# Plot all MC draws right from the start (Fig. 2.15)
# Repeat traceplots only for post-burnin draws (Fig. 2.16)
par(mfrow = c(2, 1))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we can specify the amount of 'burn-in', or how much of the chain
# we wish to discard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alpha <- out.mcmc[,1]
beta <- out.mcmc[,2]


samp <- 20
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(alpha[20] ~ samp, col = 'red', cex = 2)
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(beta[20] ~ samp, col = 'red', cex = 2)
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha[samp] + beta[samp] * x) ~ year, lwd = 3)
log.posterior(alpha[samp],beta[samp],y,x)


samp <- 200
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(alpha[samp] ~ samp, col = 'red', cex = 2)
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(beta[samp] ~ samp, col = 'red', cex = 2)
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha[samp] + beta[samp] * x) ~ year, lwd = 3)
log.posterior(alpha[samp],beta[samp],y,x)

samp <- 300
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(alpha[samp] ~ samp, col = 'red', cex = 2)
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(beta[samp] ~ samp, col = 'red', cex = 2)
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha[samp] + beta[samp] * x) ~ year, lwd = 3)
log.posterior(alpha[samp],beta[samp],y,x)


samp <- 500
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(alpha[samp] ~ samp, col = 'red', cex = 2)
plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')
points(beta[samp] ~ samp, col = 'red', cex = 2)
plot(y ~ year, 
     xlab = 'Year', ylab = 'Bee-eater counts (y)',
     pch = 21, bg = bee.cols(length(y)),
     las = 1, cex.lab = 2, type = 'b', cex = 2)
lines(exp(alpha[samp] + beta[samp] * x) ~ year, lwd = 3)
log.posterior(alpha[samp],beta[samp],y,x)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract log.posterior for every sampled pair
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lp <- NULL
for (i in 1:niter){
  lp[i] <- log.posterior(alpha[i],beta[i],y,x)
}

par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(1:niter, out.mcmc[1:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(1:niter, out.mcmc[1:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(lp, ylab = '',
     main = TeX("l($\\alpha, \\beta; y_{i}) + l(\\alpha) + l(\\beta)$"), cex.main = 2,
     xlab = "Iteration", type = 'l')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# burn the first 1000 and look at chains and log posterior
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
burn <- 1000
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(burn:niter, out.mcmc[burn:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(burn:niter, out.mcmc[burn:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(burn:niter, lp[burn:niter], ylab = '',
     main = TeX("l($\\alpha, \\beta; y_{i}) + l(\\alpha) + l(\\beta)$"), cex.main = 2,
     xlab = "Iteration", type = 'l')

par(mfrow = c(1,3), mar = c(5,5,2,2))
hist(out.mcmc[burn:niter,1], main = TeX("$\\alpha$"), breaks = 100)
hist(out.mcmc[burn:niter,2], main = TeX("$\\beta$"), breaks = 100)
hist(lp[burn:niter], main = TeX("l($\\alpha, \\beta; y_{i}) + l(\\alpha) + l(\\beta)$"), breaks = 100)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# look at individual 50 iteration runs within longer chains to think 
# about acceptance rate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(alpha[1:50] ~ c(1:50), main = TeX("$\\alpha$"), xlab = '', ylab = '')
plot(alpha[1001:1050] ~ c(1001:1050), main = TeX("$\\alpha$"), xlab = '', ylab = '')
plot(alpha[2001:2050] ~ c(2001:2050), main = TeX("$\\alpha$"), xlab = '', ylab = '')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get posterior mean, sd and 95% CRI (all post-burnin)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
post.mn <- apply(out.mcmc[burn:niter,], 2, mean) # posterior mean
post.sd <- apply(out.mcmc[burn:niter,], 2, sd) # posterior sd
CRI <- apply(out.mcmc[burn:niter,], 2,
             function(x) quantile(x, c(0.025, 0.5, 0.975))) # 95% CRI




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Acceptance rate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
acc.ratio.pb <- apply(acc[burn:niter,], 2, mean)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tuning parameters, let's change them!
# What if they're too big?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
niter <- 10000
out.mcmc <- matrix(NA, niter, 2, dimnames = list(NULL, c("alpha", "beta")))
acc <- matrix(0, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

alpha <- -1
beta <- -1 
sd.tune <- c(1, 0.1)
logpost.curr <- log.posterior(alpha, beta, y, x)

for(t in 1:niter){
  if(t %% 1000 == 0) # counter
    cat("iter", t, "\n")
  # First, update log-linear intercept (alpha)
  # ––––––––––––––––––––––––––––––––––––––––––
  # Propose candidate value of alpha
  alpha.cand <- rnorm(1, alpha, sd.tune[1]) # note tuning 'parameter'
  # Evaluate log(posterior) for proposed new alpha and
  # for current beta for the data y and covs x
  logpost.cand <- log.posterior(alpha.cand, beta, y, x)
  # Compute Metropolis acceptance ratio r
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate alpha if it meets criterion (u < r)
  if(runif(1) < r){
    alpha <- alpha.cand
    logpost.curr <- logpost.cand
    acc[t,1] <- 1 # Indicator for whether candidate alpha accepted
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Second, update log-linear slope (beta)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # propose a candidate beta-value
  beta.cand <- rnorm(1, beta, sd.tune[2])
  
  # Evaluate the log(posterior) for proposed new beta and
  # for the current alpha for the data y with covs x
  logpost.cand <- log.posterior(alpha, beta.cand, y, x)
  # Compute Metropolis acceptance ratio
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta <- beta.cand
    logpost.curr <- logpost.cand
    acc[t,2] <- 1 # Indicator for whether candidate beta accepted
  }
  out.mcmc[t,] <- c(alpha, beta) # Save samples for iteration i
  # NOTE: if proposed new values not accepted, then in this step
  # we will copy their 'old' values and insert them into object 'out'
}

lp <- NULL
for (i in 1:niter){
  lp[i] <- log.posterior(out.mcmc[i,1],out.mcmc[i,2],y,x)
}

burn <- 1000
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(burn:niter, out.mcmc[burn:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(burn:niter, out.mcmc[burn:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(burn:niter, lp[burn:niter], ylab = '',
     main = TeX("l($\\alpha, \\beta; y_{i}) + l(\\alpha) + l(\\beta)$"), cex.main = 2,
     xlab = "Iteration", type = 'l')

par(mfrow = c(1,3), mar = c(5,5,2,2))
hist(out.mcmc[burn:niter,1], main = TeX("$\\alpha$"), breaks = 100)
hist(out.mcmc[burn:niter,2], main = TeX("$\\beta$"), breaks = 100)
hist(lp[burn:niter], main = TeX("l($\\alpha, \\beta; y_{i}) + l(\\alpha) + l(\\beta)$"), breaks = 100)

acc.ratio.pb <- apply(acc[burn:niter,], 2, mean)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tuning parameters, let's change them!
# What if they're too small?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
niter <- 10000
out.mcmc <- matrix(NA, niter, 2, dimnames = list(NULL, c("alpha", "beta")))
acc <- matrix(0, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

alpha <- -1
beta <- -1 
sd.tune <- c(0.01, 0.001)
logpost.curr <- log.posterior(alpha, beta, y, x)

for(t in 1:niter){
  if(t %% 1000 == 0) # counter
    cat("iter", t, "\n")
  # First, update log-linear intercept (alpha)
  # ––––––––––––––––––––––––––––––––––––––––––
  # Propose candidate value of alpha
  alpha.cand <- rnorm(1, alpha, sd.tune[1]) # note tuning 'parameter'
  # Evaluate log(posterior) for proposed new alpha and
  # for current beta for the data y and covs x
  logpost.cand <- log.posterior(alpha.cand, beta, y, x)
  # Compute Metropolis acceptance ratio r
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate alpha if it meets criterion (u < r)
  if(runif(1) < r){
    alpha <- alpha.cand
    logpost.curr <- logpost.cand
    acc[t,1] <- 1 # Indicator for whether candidate alpha accepted
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Second, update log-linear slope (beta)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # propose a candidate beta-value
  beta.cand <- rnorm(1, beta, sd.tune[2])
  
  # Evaluate the log(posterior) for proposed new beta and
  # for the current alpha for the data y with covs x
  logpost.cand <- log.posterior(alpha, beta.cand, y, x)
  # Compute Metropolis acceptance ratio
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta <- beta.cand
    logpost.curr <- logpost.cand
    acc[t,2] <- 1 # Indicator for whether candidate beta accepted
  }
  out.mcmc[t,] <- c(alpha, beta) # Save samples for iteration i
  # NOTE: if proposed new values not accepted, then in this step
  # we will copy their 'old' values and insert them into object 'out'
}

lp <- NULL
for (i in 1:niter){
  lp[i] <- log.posterior(out.mcmc[i,1],out.mcmc[i,2],y,x)
}

burn <- 1000
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(burn:niter, out.mcmc[burn:niter,1], main = TeX("$\\alpha$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(burn:niter, out.mcmc[burn:niter,2], main = TeX("$\\beta$"), cex.main = 2,
     xlab = "Iteration", ylab = 'Posterior draw', type = 'l')

plot(burn:niter, lp[burn:niter], ylab = '',
     main = TeX("l($\\alpha, \\beta; y_{i}) + l(\\alpha) + l(\\beta)$"), cex.main = 2,
     xlab = "Iteration", type = 'l')


acc.ratio.pb <- apply(acc[burn:niter,], 2, mean)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we specify the same model in Nimble (mNM)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# And here in max-likelihood (mML)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(mML <- glm(y ~ x, family = 'poisson'))
















# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# P(>= two 6's from 3 rolled dice)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 10000
out <- NULL
for (i in 1:n){
  rolls <- rmultinom(n, 3, rep(1/6,6))
  out[i] <- length(which(rolls[6,] >= 2))/n
}

# MC approximation (works pretty well!)
mean(out)
# actual deterministic solution
true <- (1/6)^3 + 3 * ((1/6)^2 * 5/6)

par(mfrow = c(1,1), mar = c(5,5,2,2))
hist(out, breaks = 100, main = NULL, cex.lab = 1.5,
     xlab = 'Proportion of rolls with two or more sixes')
abline(v = 0.07407407, lty = 1, lwd = 3, col = 'red')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
