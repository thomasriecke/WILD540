# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we'll download packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(geoR)
library(terra)
library(fields)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set.seed: this makes everybody's 'random' numbers the same
# inv.logit: this is a function to use the logit link (we'll talk about this later)
# plotting settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)
inv.logit=plogis # redefine a function


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up an empty raster and pull coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid.res <- 50
ncell <- grid.res^2
r=rast(ncol=grid.res,nrow=grid.res,extent=c(-grid.res,grid.res,-grid.res,grid.res))                      # make an empty raster
r[]=0                                                                # cero; 0
s.r=crds(r)   

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create two random fields
# one will be quite homogenous, we'll call this zS (smooth)
# one will be quite variable, we'll call this zV (variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z1 <- r
z2 <- r
z1[] <- grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(0.1,10), mean = 0)$data 
z2[] <- grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10), mean = 0)$data 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we'll make a (somewhat hasty) figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
fields::image.plot(z1, zlim = c(-5,5), main = 'Fairly Homogenous', axes = 0)
fields::image.plot(z2, zlim = c(-5,5), main = 'Heterogenous and clumped', axes = 0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's look at the true means and variances for each landscape
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean(z1[])
sd(z1[])
mean(z2[])
sd(z2[])




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's sample the two sites (20 cells each)
# and see if they're 'statistically different'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n1 <- 20
n2 <- 20

s1 <- sample(ncell, n1, replace = F)
s2 <- sample(ncell, n2, replace = F)

y1 <- z1[][s1]
y2 <- z2[][s2]

test <- t.test(y1,y2)
str(test)
test

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's do the same thing 10,000 times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nSims <- 10000
n1 <- 20
n2 <- 20
p <- rep(NA, nSims)
for (ii in 1:nSims){
  s1 <- sample(ncell, n1, replace = F)
  s2 <- sample(ncell, n2, replace = F)
  
  y1 <- z1[][s1]
  y2 <- z2[][s2]
  
  test <- t.test(y1,y2)
  p[ii] <- test$p.value
}
par(mfrow = c(1,1))
hist(p, xlab = 'p-values', main = NULL, breaks = 100)
length(which(p < 0.05))/nSims
length(which(p < 0.01))/nSims

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's do the same thing 10,000 times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nSims <- 10000
n1 <- 50
n2 <- 50
p <- rep(NA, nSims)
for (ii in 1:nSims){
  s1 <- sample(ncell, n1, replace = F)
  s2 <- sample(ncell, n2, replace = F)
  
  y1 <- z1[][s1]
  y2 <- z2[][s2]
  
  test <- t.test(y1,y2)
  p[ii] <- test$p.value
}
par(mfrow = c(1,1))
hist(p, xlab = 'p-values', main = NULL, breaks = 100)
length(which(p < 0.05))/nSims
length(which(p < 0.01))/nSims

