# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(terra)
library(geoR)
library(latex2exp)

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate values (v) to fill grid (truth is z)
# more fun Gaussian random field reading (not required)
# https://www.paulamoraga.com/book-spatial/geostatistical-data-1.html
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
v <- grf(1, grid = st_coordinates(p), cov.pars=c(1,5), mean = 0)$data
z <- st_as_sf(g)
z$v <- v

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's sample (s) the grid given a sample size... 
# n: number of samples
# s: the actual samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first, we'll make a plot of our landscape
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mar = c(0,0,0,0), # defines inner margins of plot
    oma = c(1,1,1,5)) # defines outer margins of plot
plot(z, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# then we'll sample it and add our random points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 25                            # sample size
s <- sample(nc, n, replace = F)    # randomly sample n cells
s                                  # these are the randomly selected n grid cells
points(p[s], col = adjustcolor('black', alpha = 0.5), pch = 19, 
       xlim = st_bbox(z), ylim = st_bbox(z))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summary statistics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean(v)     # mean of the entire landscape
mean(v[s])  # mean of the sample
sd(v)       # sd of the entire landscape
sd(v[s])    # sd of the sample
sd(v[s])/sqrt(n)    # se


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# so... we could keep doing that by hand 1000 times... or we could put
# it in a loop.... so... loop?!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# number of simulated samples
nSims <- 1000

# these are empty vectors to store data (one value for each iteration)
mu <- rep(NA, nSims)
stdev <- rep(NA, nSims)
se <- rep(NA, nSims)
out <- rep(NA, nSims)
# sample size
n <- 25                            

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# and here's a for-loop
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in 1:nSims){
  # randomly sample n cells
  s <- sample(nc, n, replace = F)    
  # mean of the sample
  mu[i] <- mean(v[s]) 
  # standard deviation of the sample
  stdev[i] <- sd(v[s]) 
  # let's save the standard error of the mean as well
  se[i] <- stdev[i]/sqrt(n)
  # is the true mean outside of the confidence intervals?
  out[i] <- (mu[i]+1.96*se[i]) < mean(v) | (mu[i]-1.96*se[i]) > mean(v)
}

par(mfrow = c(1,2), mar = c(5,5,2,2), oma = rep(0,4))
hist(mu, main = NULL, xlab = TeX("Sample means ($\\mu$)"), breaks = 50)
abline(v = mean(v), col = 'red', lty = 1, lwd = 3)
hist(stdev, main = NULL, xlab = TeX("Sample st. dev.'s ($\\sigma$)"), breaks = 50)
abline(v = sd(v), col = 'red', lty = 1, lwd = 3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 1
# how many samples and 95% CIs contain the true mean?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(out)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This was exactly 950 when I ran it (out of 1000) 
# it is (on average) a little less than that with this simulation
# why?
# What are 95% confidence intervals?
# "Were this procedure to be repeated on numerous samples, the proportion of calculated 95% 
# confidence intervals that encompassed the true value of the population parameter would tend toward 95%."[Cox 1974]
#
# Why might we be a little less than 95%?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 2
# What if we're interested in testing different sample sizes?
# We could re-run the code above with a different sample size (n)
# or we could modify the code to have n sample sizes
# talk to your neighbor, then review the code below
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# number of simulated samples
nSims <- 1000
# sample size
n <- c(25,50)  

# these are empty matrices to store data (one value for each iteration)
mu <- matrix(NA, nSims, length(n))
stdev <- matrix(NA, nSims, length(n))
se <- matrix(NA, nSims, length(n))
out <- matrix(NA, nSims, length(n))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# and here's a for-loop
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in 1:nSims){
  for (j in 1:length(n)){
    # randomly sample n cells
    s <- sample(nc, n[j], replace = F)    
    # mean of the sample
    mu[i,j] <- mean(v[s]) 
    # standard deviation of the sample
    stdev[i,j] <- sd(v[s]) 
    # let's save the standard error of the mean as well
    se[i,j] <- stdev[i]/sqrt(n[j])
    
    out[i,j] <- (mu[i,j]+1.96*se[i,j]) < mean(v) | (mu[i,j]-1.96*se[i,j]) > mean(v)
  }
}

par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
boxplot(mu, names = n, xlab = 'Sample size', ylab = TeX("Mean ($\\mu$)"),
        cex.lab = 2, las = 1)
boxplot(stdev, names = n, xlab = 'Sample size', ylab = TeX("St. dev. ($\\sigma$)"),
        cex.lab = 2, las = 1)
boxplot(se, names = n, xlab = 'Sample size', ylab = TeX("Standard error"),
        cex.lab = 2, las = 1)

colMeans(out)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 3
# Add another sample size (100?) to the code above and remake the figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 4 [see slide 7]
# How would we determine how many of our samples were within x
# of the true mean? 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prop <- NULL
diff <- 0.2
# could this be in a loop?
prop[1] <- length(which(abs(mu[,1] - mean(v)) < diff))/nSims
prop[2] <- length(which(abs(mu[,2] - mean(v)) < diff))/nSims

barplot(prop, ylab = 'Proportion of samples within 0.1 of truth',
        las = 1, cex.lab = 1.5, ylim = c(0,1),
        xlab = 'Sample size', names = n)
box()
