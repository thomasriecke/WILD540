# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Class script 2:
# an introduction to power analysis and random sampling...
# We'll simulate two landscapes
# then we'll sample them (many times)
#
# first, we'll download some packages
# the spatial aspects will clean up as we go...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(geoR)
library(terra)
library(fields)
library(stars)
library(latex2exp)
library(sf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set.seed: this makes everybody's 'random' numbers the same
# you don't have to run this, but it can be helpful at times, particularly
# if you're archiving code...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set.seed(123)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up an empty SpatRaster and pull coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid.res <- 50
ncell <- grid.res^2
r=rast(ncol=grid.res,nrow=grid.res,extent=c(-grid.res,grid.res,-grid.res,grid.res))                      # make an empty raster
r[]=0                                                                # cero; 0
sr=crds(r)   

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ok... at this point we have four things in our environment
# grid.res: the dimensions (one-side) of each grid = 50
# ncell: the total number of cells in the grid = 50 * 50 = 2500
# r: an empty grid that is 50 * 50 where every cell is 0
# sr: the spatial locations of each cell in the grid (like coordinates or UTMs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll create two random landscapes
# one will be quite homogenous, we'll call this z1 (smooth)
# one will be quite variable, we'll call this z2 (variable)
# here is a really nice write-up on Gaussian random fields (grf):
# https://structures.uni-heidelberg.de/blog/posts/gaussian-random-fields/index.php
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z1 <- r
z2 <- r
z1[] <- grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(0.1,10), mean = 0)$data 
z2[] <- grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10), mean = 0)$data 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we'll make a (somewhat hasty) figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(2.1,2.1,2.1,2.1))
fields::image.plot(z1, zlim = c(-5,5), main = 'Fairly Homogenous', axes = 0)
fields::image.plot(z2, zlim = c(-5,5), main = 'Heterogenous and clumped', axes = 0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's look at the true means and 
# variance (actually st. dev's) for each landscape
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean(z1[])
sd(z1[])
mean(z2[])
sd(z2[])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's sample the two sites (20 cells each)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n1 <- 20
n2 <- 20
s1 <- sample(ncell, n1, replace = F)
s2 <- sample(ncell, n2, replace = F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# these are the 20 randomly selected points from each site
# note that they're just numbered according to the order 
# of the grid (so 1 would be upper left cell, 2500 would be lower right cell)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s1
s2


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fields::image.plot(z1, zlim = c(-5,5), main = 'Fairly Homogenous', axes = 0)
points(sr[s1,2] ~ sr[s1,1])

fields::image.plot(z2, zlim = c(-5,5), main = 'Heterogenous and clumped', axes = 0)
points(sr[s2,2] ~ sr[s2,1])

y1 <- z1[][s1]
y2 <- z2[][s2]



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# introducing for-loops
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in 1:10){
  print(i)
}

for (i in 1:10){
  print(i + 3)
}

for (i in 1:10){
  print(i * 3)
}

text.string <- c('I','am','not','entirely','convinced',
                 'regarding','the','utility','of','for-loops')
for (i in 1:10){
  print(text.string[i])
}

# We can use vectors to index
loop <- 1:10
for (i in loop){
  print(text.string[i])
}

# we can have nested for-loops
loop <- 1:10
for (i in loop){
  print(text.string[i])
  for (j in 1:3){
    print(j)
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# So...
# for-loops do everything inside of them however many times
# we tell them to
# 
# We can save that output as we go. That can be incredibly useful!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's use a for-loop to sample 
# our landscape 1000 times...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nSims <- 1000

# we'll pull 20 points each time...
n1 <- 20
n2 <- 20

# create vectors to store means and standard deviations
m1 <- rep(NA, nSims)
m2 <- rep(NA, nSims)
# 
s1 <- rep(NA, nSims)
s2 <- rep(NA, nSims)

for (ii in 1:nSims){
  samp1 <- sample(ncell, n1, replace = F)
  samp2 <- sample(ncell, n2, replace = F)
  
  y1 <- z1[][samp1]
  y2 <- z2[][samp2]
  
  m1[ii] <- mean(y1)
  m2[ii] <- mean(y2)
  
  s1[ii] <- sd(y1)
  s2[ii] <- y2
  
}
par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
hist(m1, xlab = TeX("Estimates of mean ($\\mu_1$)"), 
     main = NULL, breaks = 100, xlim = c(-1,1))
hist(m2, xlab = TeX("Estimates of mean ($\\mu_2$)"), 
     main = NULL, breaks = 100, xlim = c(-1,1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 1
# what if our goal was for our estimate of the mean to be within
# 0.25 of truth 90% of the time?
#
# are we there yet for landscape 1?
# what about for landscape 2?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's use a for-loop to sample 
# our landscape 1000 times with different sample sizes?!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nSims <- 1000

# we'll pull 20 points each time...
nBlocks <- 2
n1 <- c(20, 50)
n2 <- c(20, 50)

# create matrices to store means and standard deviations
m1 <- matrix(NA, nSims, nBlocks)
m2 <- matrix(NA, nSims, nBlocks)
# 
s1 <- matrix(NA, nSims, nBlocks)
s2 <- matrix(NA, nSims, nBlocks)

for (ii in 1:nSims){
  for (jj in 1:nBlocks){
    samp1 <- sample(ncell, n1[jj], replace = F)
    samp2 <- sample(ncell, n2[jj], replace = F)
    
    y1 <- z1[][samp1]
    y2 <- z2[][samp2]
    
    m1[ii,jj] <- mean(y1)
    m2[ii,jj] <- mean(y2)
    
    s1[ii,jj] <- sd(y1)
    s2[ii,jj] <- sd(y2)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
boxplot(m1,  
     main = NULL, breaks = 100,
     ylab = TeX("Estimates of mean ($\\mu_1$)"), xlab = 'Sample size',
     names = n1)
boxplot(m2,  
        main = NULL, breaks = 100,
        ylab = TeX("Estimates of mean ($\\mu_2$)"), xlab = 'Sample size',
        names = n1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 2
# What if we wanted to know what proportion of our estimates for each
# sample size were within 0.1 of truth?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM (really just think about this!) 4
# Here we'll make a box within our grid, i.e., we'll drop a 'study area'
# into our landscape. Play with the dimensions a bit, changing the
# shape of the polygong
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sr[1:10,] # note that the grid coordinates are for cell centers

# here we'll make a grid
sa.grid <- 30
min <- -sa.grid
max <- sa.grid
study <- list(matrix(
  c(min,min, max,min, max,max, min,max, min,min), 
  ncol = 2, byrow = T))
study <- sf::st_polygon(study)

# let's plot it.. cool!
fields::image.plot(z1, zlim = c(-5,5), main = 'Fairly Homogenous', axes = 0)
plot(study, add = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# there are cool ways to do this using spatial R packages
# but this is a nice and simple way with logical statements
# 
# in the future we'll create multiple study areas within landscapes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inside <- which(sr[,1] > -sa.grid & sr[,1] < sa.grid &
                sr[,2] > -sa.grid & sr[,2] < sa.grid)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we could subsample this smaller grid (inside)
# within our broader landscapes (z1 and z2)
#
# we'll look at different ways to do this next week...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



