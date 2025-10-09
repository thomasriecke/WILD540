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


# random spatiotemporal variation
epsilon <- MASS::mvrnorm(1, rep(0, nc), Sigma)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# randomly sample 100 sites prior to impact
# resample the same sites after impact occurs across half of study area
# sample 25 for treatments (e.g., exclosures)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n  <- 100 
# total sample size (we would nest equal sample sizes, but we'll use)
s <- sample(nc, n, replace = F) # s is our sample

# assign values to sampled cells
sites <- rep(s, 2)
aft <- c(rep(0,n), rep(1,n))
imp <- rep(0, n*2)
imp[which(sites %in% blk)] <- 1

data <- data.frame(cbind(sites, aft, imp))





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign values (v) to grids before (z1) and after (z2) treatment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z1 <- st_as_sf(g)
z2 <- st_as_sf(g)


# here we generate two empty vectors of length 2500 (a value for each cell)
# we'll fill some of these with 1's if they're 'in' or 'treated'
trt <- rep(0, nc)
ar <- rep(0, nc)
ar[blk] <- 1; table(ar)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 1a: here we changed effects
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# effect of treatment (beta[1]) and year (beta[2])
beta <- c(0.5, 0.5)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 3: use an exp() function to only have positive values
# density?
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# linear model
z1$v <- exp(epsilon)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 1b: here we added a little variance between before and after
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
z2$v <- exp(epsilon + beta[1] * ar + beta[2] + rnorm(nc, 0, 0.1))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make figures
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 3: note we changed the breaks to only be > 0
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(z1, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,25,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
# let's add our sampled and treated points (treated points will be dark blue)
points(p[s], cex = 0.5, pch = 19)


plot(z2, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,25,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
# let's add our sampled and treated points (treated points will be dark blue)
points(p[s], cex = 0.5, pch = 19)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract data from our sampled plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data$y <- NA
str(data)
data$y[aft == 0] <- z1$v[s]
data$y[aft == 1] <- z2$v[s]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new, we'll make a boxplot plot of our estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5,5,2,2), oma = c(0,0,0,0))
boxplot(data$y ~ data$imp + data$aft , ylab = 'Data',
        names = c('Before:Control', 'Before:Impact', 'After:Control', 'After:Impact'),
        xlab = 'Group', col = c('white','red','white','red'), cex.lab = 2,
        las = 1)







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 3: now we model the log
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
str(data)
m <- lm(log(y) ~ imp + aft + imp:aft, data = data)
summary(m)

# can also be written as 
m <- lm(log(y) ~ imp*aft, data = data)
summary(m)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# inference
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# average for control prior to impact was -0.8115
# no significant difference between control and impact areas prior to impact
# effect of year change (after) was 1
# difference in change between control group and impact group was 1...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# obtain predictions from the model and make a plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 3: we plot exponents
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# obtain predictions from 'new data'
newdata <- data.frame(aft = c(0,0,1,1), imp = c(0,1,0,1))
p <- predict(m, newdata = newdata, interval = 'confidence')
# what's going on here?
p <- data.frame(p)
plot(exp(p$fit), xaxt = 'n', ylim = c(0,6), xlim = c(0.5,4.5),
     xlab = 'Group', ylab = 'Estimate', cex.lab = 2)
axis(side = 1, at = 1:4, labels = c('control','impact','control','impact'))

# add actual data points
points(data$y[data$imp == 0 & data$aft == 0] ~ 
         jitter(rep(1,length(data$y[data$imp == 0 & data$aft == 0]))), cex = 0.5)
points(data$y[data$imp == 1 & data$aft == 0] ~ 
         jitter(rep(2,length(data$y[data$imp == 1 & data$aft == 0]))), cex = 0.5)
points(data$y[data$imp == 0 & data$aft == 1] ~ 
         jitter(rep(3,length(data$y[data$imp == 0 & data$aft == 1]))), cex = 0.5)
points(data$y[data$imp == 1 & data$aft == 1] ~ 
         jitter(rep(4,length(data$y[data$imp == 1 & data$aft == 1]))), cex = 0.5)

# add confidence intervals
arrows(1:4, exp(p$lwr), 1:4, exp(p$upr), length = 0, lwd = 4)
# add means
points(exp(p$fit), pch = 21, bg = c('white','red','white','red'), cex = 2.5)
abline(v = 2.5, lty = 2)
text(1.5,6,'Before')
text(3.5,6,'After')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  # Priors 
  alpha ~ dnorm(0, sd = 10)
  beta[1] ~ dnorm(0, sd = 10)
  beta[2] ~ dnorm(0, sd = 10)
  beta[3] ~ dnorm(0, sd = 10)  
  sigma ~ dgamma(1,1)
  # Likelihood
  for(i in 1:n){
    # mu: expected value
    mu[i] <- alpha + beta[1] * imp[i] + beta[2] * aft[i] + beta[3] * imp[i] * aft[i]
    # data are a function of expected value and error
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
} )

# Bundle data (same as for JAGS)
str(dataList <- list(y = log(data$y), n = nrow(data), 
                     imp = data$imp, aft = data$aft)) # not shown


# Can use same function to generate starting values as for JAGS
inits <- function()
  list (alpha = 0, beta = c(0,0,0), sigma = 2)

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's make a plot from our nimble output
# first we'll put all the chains together into one big chain
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- do.call(rbind.data.frame, m$samples)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# then we'll create a matrix to store predictions
# 4 columns, one for each group (out, control; out, trt; in, control; in, trt)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pr <- matrix(NA, nrow(tmp), 4)
q <- matrix(NA, 5, 4)
pr[,1] <- tmp$alpha
pr[,2] <- tmp$alpha + tmp$`beta[1]`
pr[,3] <- tmp$alpha + tmp$`beta[2]`
pr[,4] <- tmp$alpha + tmp$`beta[1]` + tmp$`beta[2]` + tmp$'beta[3]'
for (j in 1:4){
  q[,j] <- quantile(pr[,j], c(0.025,0.25,0.5,0.75,0.975))
}

par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(NA, xaxt = 'n', ylim = c(0,6), xlim = c(0.5,4.5),
     xlab = 'Group', ylab = 'Estimate', cex.lab = 2)
axis(side = 1, at = 1:4, labels = c('control','Impact','control','Impact'))

# add actual data points
points(data$y[data$imp == 0 & data$aft == 0] ~ 
         jitter(rep(1,length(data$y[data$imp == 0 & data$aft == 0]))), cex = 0.5)
points(data$y[data$imp == 1 & data$aft == 0] ~ 
         jitter(rep(2,length(data$y[data$imp == 1 & data$aft == 0]))), cex = 0.5)
points(data$y[data$imp == 0 & data$aft == 1] ~ 
         jitter(rep(3,length(data$y[data$imp == 0 & data$aft == 1]))), cex = 0.5)
points(data$y[data$imp == 1 & data$aft == 1] ~ 
         jitter(rep(4,length(data$y[data$imp == 1 & data$aft == 1]))), cex = 0.5)

# add confidence intervals
arrows(1:4, exp(q[1,]), 1:4, exp(q[5,]), length = 0, lwd = 4)
# add means
points(exp(q[3,]), pch = 21, bg = c('white','red','white','red'), cex = 2.5)
abline(v = 2.5, lty = 2)
text(1.5,6,'Before')
text(3.5,6,'After')









# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# additional plotting code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(NA, xaxt = 'n', ylim = c(0,6), xlim = c(0.5,2.5),
     xlab = 'Time', ylab = 'Estimate', cex.lab = 2)
axis(side = 1, at = 1:2, labels = c('Before','After'))

# add actual data points
points(data$y[data$imp == 0 & data$aft == 0] ~ 
         jitter(rep(0.9,length(data$y[data$imp == 0 & data$aft == 0]))), cex = 0.5)
points(data$y[data$imp == 1 & data$aft == 0] ~ 
         jitter(rep(1.1,length(data$y[data$imp == 1 & data$aft == 0]))), cex = 0.5)
points(data$y[data$imp == 0 & data$aft == 1] ~ 
         jitter(rep(1.9,length(data$y[data$imp == 0 & data$aft == 1]))), cex = 0.5)
points(data$y[data$imp == 1 & data$aft == 1] ~ 
         jitter(rep(2.1,length(data$y[data$imp == 1 & data$aft == 1]))), cex = 0.5)

# add confidence intervals
arrows(c(0.9,1.1,1.9,2.1), exp(q[1,]), c(0.9,1.1,1.9,2.1), exp(q[5,]), length = 0, lwd = 4)
# add means
points(exp(q[3,]) ~ c(0.9,1.1,1.9,2.1), pch = 21, bg = c('white','red','white','red'), cex = 2.5)
legend('topleft', legend = c('Control','Impact'), box.lty = 0,
       pch = c(21,21), pt.bg = c('white','red'), cex = 1.5)
box()





