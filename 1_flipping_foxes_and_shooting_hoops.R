# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMPLE 1: Sampling from a binomial distribution
# first we'll define the probability of survival and call it 'phi'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phi <- 0.5 # probability of survival
phi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we're going to use a couple of functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
?rbinom()
?pbinom()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll see if one fox survives (1: alive, 0: not alive)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rbinom(1, 1, 0.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll see how many out of 100 foxes survive
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rbinom(1, 100, phi)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now see how many of 100 foxes (f) at 100 sites (n) survive
# this would be the most ambitious study of foxes in history
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100  # number of sites
f <- 100  # number of foxes
d <- rbinom(n, f, phi)
d

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's stop a moment and think, it would be easy to look at the lowest
# and highest number of survivors and think that these are particularly
# poor and good sites, respectively.
#
# but we know that's not true! we just made these data up, and each site
# had the exact same survival probability
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we'll calculate the actual expected mean and variance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# observed mean
sum(d)/n
mean(d)
# expected mena
phi * n


# observed variance
sum((d - mean(d))^2)/(n-1)
var(d)
sd(d)

# expected variance
n * phi * (1-phi)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# play with the number of trials (n) and the probability (p)
# below and see what happens visually (via the histogram),
# to the mean, and to the standard deviation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100
p <- 0.1
hist(rbinom(10000, n, p), breaks = 50,
     main = paste0('Binomial(',n,', ',p,')'),
     xlab = '')

# print mean
n * p
# print standard deviation
sqrt(n * p * (1-p))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's think about possible potential outcomes!
# we can calculate these by hand... here is the equation for the 
# probability of 50 foxes surviving
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(factorial(n)/(factorial(100-50) * factorial(50))) * (0.5^50) * (1 - 0.5)^(100-50)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# work with your neighbor to adjust the above equation to calculate the 
# probability of 49 foxes surviving
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Want to do this for all the numbers instead of one at a time?
# dbinom will give us the probability of each possible outcome
# anywhere from 0 to 100 foxes surviving
# we'll call this 'pmf' for probability mass function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pmf <- dbinom(0:100, f, phi)
pmf[51] # same as above
plot(pmf, type = 'l', xlab = 'Number of surviving foxes given n = 100 and p = 0.5', ylab = 'P(Outcome)')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# explore the properties of the Poisson distribution similar to 
# problem 1 above
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hist(rpois(10000, 4), main = NULL, xlab = 'Poisson(4)')
?rpois()
?dpois()





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMPLE 2: the hot hand fallacy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sims <- 10000
result <- NULL
for (ii in 1:sims){
  phi <- 0.5 # probability of making a shot
  
  # take a bunch of shots (n = 100)
  n <- 100
  d <- rbinom(n, 1, phi)
  
  # let's see what the average shooting percentage was
  mean(d)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # now we'll find all the shots that occurred after three makes (u)
  # we'll use a bit of perl (gregexpr; perl is a different coding language
  # that can also be used in R) here that looks for strings of '111's in
  # the data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp <- paste0(d[1:(n-1)], collapse = '')
  h <- gregexpr('1(?=11)', tmp, perl = T)[[1]] + 3
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # h is now a vector of all the shots that come after three 1's in a row
  # now we'll see if the average of the shots after the shooter is 
  # hot (3 makes in a row; d[h]),  is worse than the total mean (d)...
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  result[ii] <- mean(d) > mean(d[h])
}
table(result)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the average shooting percentage (mean(d)) is greater than the 
# shooting percentage during a 'hot' streak in 58% of simulations,
# even though
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~











