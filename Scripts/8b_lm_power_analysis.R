# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# There are two power analyses in this script:
#
# In the first, we randomly generate beta values (effect sizes)
# with a mean of 1 and a standard deviation of 0.5
#
# In the second, we supply fixed effect sizes
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# POWER ANALYSIS 1
# This is a power analysis given randomly generated covariate values (x),
# randomly generated effect sizes (beta), 
# and data that are simulated as a function of the simulated covariate
# values, simulated effect sizes (beta), and user-defined sample sizes (n)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# number of simulations
nSims <- 1000

# sample size
n <- c(5,50,500)

# residual variance
sigma <- 1

# arrays to store b and p
# beta: simulated effect size
beta <- array(NA, dim = c(nSims, length(n)))

# b: estimate of the beta value, or effect size, from a linear model, lm()
b <- array(NA, dim = c(nSims, length(n)))
# p: p-value for b
p <- array(NA, dim = c(nSims, length(n)))

for (i in 1:nSims){
  # print(i)
  for (j in 1:length(n)){
      # new sample size (we can adjust based on stochasticity or failure rate)
      new <- n[j]
      
      # random effect size
      beta[i,j] <- rnorm(1, 1, 0.5)
      
      # x; covariate value, randomly generated
      x <- rnorm(new, 0, 1)
      
      # y: data, created as a function of x
      y <- rnorm(new, beta[i,j] * x, sigma)
      
      # this fits the model
      # y ~ normal(intercept + slope * x, sigma)
      # we'll store the estimate of the slope (b)
      # and the p-value for the slope (p)
      m <- lm(y ~ x)
      b[i,j] <- m$coefficients[2]
      p[i,j] <- summary(m)$coefficients[2,4]
      
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a grid plot of beta estimates given different effect and sample sizes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,length(n)), oma = c(2,2,0,0))
for (k in 1:length(n)){
  plot(b[,k] ~ beta[,k], names = n, main = paste0("n", '=', n[k]))
}
mtext('Estimated beta', side = 2, outer = T)
mtext('Simulated beta', side = 1, outer = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a grid plot of p-values given different effect and sample sizes
# why the bump around 0?!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,length(n)))
for (k in 1:length(n)){
  plot(p[,k] ~ beta[,k], names = n, main = paste0("n", '=', n[k]),
       ylim = c(0,1))
}
mtext('p-value', side = 2, outer = T)
mtext('Simulated absolute value of beta', side = 1, outer = T)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate 'power'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pow <- rep(NA, length(j))
alpha <- 0.05 # significance level
# for each sample size (j) and effect size (k)
for (j in 1:length(n)){
    # how many simulations had a 1) positive effect
    #                            2) p-value less than 0.05
    #                            3) beta value within 0.2 of the true value
    pow[j] = length(which(p[,j] < alpha &
                            abs(b[,j] - beta[,j]) < 0.2))/nSims
}

par(mfrow = c(1,1))
barplot(pow, names = n, main = NULL, ylim = c(0,1))
mtext('Power', side = 2, outer = T)
mtext('Sample size', side = 1, outer = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# with a sample size of 500, we have 95% power to detect a significant effect
# of any effect size...
# 
# with a sample size of 50, we have about 75% power...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# POWER ANALYSIS 2
# This is a power analysis given randomly generated covariate values (x)
# and data that are simulated as a function of the simulated covariate
# values, user-defined effect sizes (beta), and user-defined sample sizes (n)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# number of simulations
nSims <- 100

# sample size
n <- c(10,30,50,100,250)

# effect size
beta <- c(0.1, 0.5, 1)

# residual variance
sigma <- 1

# arrays to store b and p
# b: estimate of the beta value, or effect size, from a linear model, lm()
b <- array(NA, dim = c(nSims, length(n), length(beta)))
# p: p-value for b
p <- array(NA, dim = c(nSims, length(n), length(beta)))

for (i in 1:nSims){
  print(i)
  for (j in 1:length(n)){
    for (k in 1:length(beta)){
      
      # new sample size (we can adjust based on stochasticity or failure rate)
      new <- n[j]
      
      # x; covariate value, randomly generated
      x <- rnorm(new, 0, 1)
      # y: data, created as a function of x
      y <- rnorm(new, beta[k] * x, sigma)
      
      # this fits the model
      # y ~ normal(intercept + slope * x, sigma)
      # we'll store the estimate of the slope (b)
      # and the p-value for the slope (p)
      m <- lm(y ~ x)
      b[i,j,k] <- m$coefficients[2]
      p[i,j,k] <- summary(m)$coefficients[2,4]
      
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a grid plot of beta estimates given different effect and sample sizes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,length(beta)), oma = c(2,2,0,0))
for (k in 1:length(beta)){
  vioplot::vioplot(b[,,k], names = n, main = paste0("beta", '=', beta[k]))
}
mtext('Beta', side = 2, outer = T)
mtext('Sample size', side = 1, outer = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a grid plot of p-values given different effect and sample sizes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,length(beta)))
for (k in 1:length(beta)){
  vioplot::vioplot(p[,,k], names = n, main = paste0("beta", '=', beta[k]))
}
mtext('p-value', side = 2, outer = T)
mtext('Sample size', side = 1, outer = T)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate 'power'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pow <- matrix(NA, length(n), length(beta))
alpha <- 0.05 # significance level
# for each sample size (j) and effect size (k)
for (j in 1:length(n)){
  for (k in 1:length(beta)){
    # how many simulations had a 1) positive effect
    #                            2) p-value less than 0.05
    #                            3) beta value within 0.2 of the true value
    pow[j,k] = length(which(b[,j,k] > 0 &
                              p[,j,k] < alpha &
                              abs(b[,j,k] - beta[k]) < 0.2))/nSims
  }
}

par(mfrow = c(1,length(beta)))
for (k in 1:length(beta)){
  barplot(pow[,k], names = n, main = paste0("beta", '=', beta[k]), ylim = c(0,1))
}
mtext('Power', side = 2, outer = T)
mtext('Sample size', side = 1, outer = T)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# When beta = 1, sample sizes of 50 and 100 had a power > 0.8
# when beta = 0.5, 50 is borderline and 100 is ok, 250 is great for both
# when beta = 0.1, a sample size of 100 only had ~20% power... a sample size
# of 250 only had ~40% power.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# So... the 'golden' sample size depends on the effect size and our decisions 
# for how close our parameter estimate needs to be to truth
# these will always change based on our study design and system
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Things we learned from this:
# -if our effect size (effect of x on y) is 0.1, we'll need a HUGE number of samples
#  to learn anything. If it's larger, we could probably truncate around 50 samples
#
# -there are diminishing returns. The difference in power b/w 100 and 250 samples is relatively
#  
# - we could incorporate other sources of 'noise' such as sampling failure to 
#   account for potential problems in the field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
