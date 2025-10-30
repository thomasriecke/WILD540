# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(terra)
library(latex2exp)
library(nimble)
library(ASMbook)
library(lme4)
library(MuMIn)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# contents
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) simulate a remotely-sensed landscape {Lines 20-100}
# (2) simulate abundance as a function of landscape (Lines 100-120)
# (3) sample abundance (Lines 100-120)
# (4) run four models and calculate AICc, think about model weights
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) simulate a landscape {Lines 20-100}
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use a geospatial model to generate something similar to a GRF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# D: CREATE A DISTANCE MATRIX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D <- st_distance(p, p)
dim(D)
par(mfrow = c(1,1), mar = c(5,5,2,2))
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use Cholesky decomposition to speed up this line of code
# epsilon <- MASS::mvrnorm(1, rep(0, nc), Sigma)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L <- chol(Sigma)
dim(Sigma)
x <- rnorm(nc, 0, 1)
soil <-  as.numeric(x %*% L)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign values to z
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- st_as_sf(g)
z$soil <- soil
plot(z['soil'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a patchy burn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
b <- rnorm(nc, 0, 1)
risk <- as.numeric(b %*% L)
burn <- rbinom(nc, 1, plogis(risk))
z$burn <- burn
plot(z['burn'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,1,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (2) simulate abundance as a function of landscape (Lines 100-130)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(1, 0.25, 0.45)
E <- exp(alpha[1] + alpha[2] * soil + alpha[3] * burn + rnorm(nc, 0, 0.1))
n <- rpois(nc, E)
y <- rpois(nc, E)
z$n <- n
z$y <- y
plot(z["y"], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(0,20,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's plot our true outcomes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(jitter(n) ~ soil, data = data, ylab = 'Abundance (n)', xlab = 'Soil richness',
     col = burn + 1, cex.lab = 2, pch = 19)

boxplot(n ~ burn, col = c('grey', 'firebrick2'), names = c('Unburned', 'Burned'),
        xlab = '', ylab = 'Abundance (n)', cex.lab = 2)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample the landscape (counts)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S <- 100
s <- sample(nc, S, replace = F)

plot(z['soil'], 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
points(p[s], pch = 19)


data <- data.frame(matrix(NA, nrow = S, ncol = 3))
names(data) <- c('y','s','b')
data$y <- z$y[s]
data$s <- z$soil[s]
data$b <- z$burn[s]

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(y ~ s, data = data, ylab = 'Count (y)', xlab = 'Soil richness',
     col = b + 1, cex.lab = 2, pch = 19)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (4) run a model set!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a) ML estimates and predictions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(m0 <- glm(y ~ 1, data = data, family = 'poisson'))
summary(m1 <- glm(y ~ s, data = data, family = 'poisson'))
summary(m2 <- glm(y ~ b, data = data, family = 'poisson'))
summary(m3 <- glm(y ~ s + b, data = data, family = 'poisson'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model selection
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AICc(m0,m1,m2,m3)
model.sel(m0,m1,m2,m3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's calculate the likelihood by hand
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p0 <- predict(m0)
p1 <- predict(m1)
p2 <- predict(m2)
p3 <- predict(m3)

# log-likelihoods
ll0 <-  -sum(exp(p0)) + sum(data$y * p0) - sum(log(factorial(data$y)))
ll1 <-  -sum(exp(p1)) + sum(data$y * p1) - sum(log(factorial(data$y)))
ll2 <-  -sum(exp(p2)) + sum(data$y * p2) - sum(log(factorial(data$y)))
ll3 <-  -sum(exp(p3)) + sum(data$y * p3) - sum(log(factorial(data$y)))


# likelihood for saturated model (lu)
lu <-  -sum(data$y[data$y > 0]) + sum(data$y[data$y > 0] * log(data$y[data$y > 0])) - sum(log(factorial(data$y[data$y > 0])))

# residual deviance
d0 <- -2 * (ll0 - lu)
d1 <- -2 * (ll1 - lu)
d2 <- -2 * (ll2 - lu)
d3 <- -2 * (ll3 - lu)

# calculate AIC
aic0 <- length(m0$coefficients) * 2 + -2*ll0
aic1 <- length(m1$coefficients) * 2 + -2*ll1 
aic2 <- length(m2$coefficients) * 2 + -2*ll2 
aic3 <- length(m3$coefficients) * 2 + -2*ll3 



# calculate AICc
n.params <- c(1,2,2,3)
AIC(m0,m1,m2,m3)
AICc(m0,m1,m2,m3)
k0 <- 1; k1 <- 2; k2 <- 2; k3 = 3
aicc0 <- aic0 + (2*k0^2 + 2*k0)/(nrow(data) - k0 - 1)
aicc1 <- aic1 + (2*k1^2 + 2*k1)/(nrow(data) - k1 - 1)
aicc2 <- aic2 + (2*k2^2 + 2*k2)/(nrow(data) - k2 - 1)
aicc3 <- aic3 + (2*k3^2 + 2*k3)/(nrow(data) - k3 - 1) 



# calculate relative likelihoods and model weights
# note these equations will only work if mod3 is the top model
rl0 <- exp(-0.5 * (aic0 - aic3))
rl1 <- exp(-0.5 * (aic1 - aic3))
rl2 <- exp(-0.5 * (aic2 - aic3))
rl3 <- exp(-0.5 * (aic3 - aic3))

w0 <- rl0/(rl0 + rl1 + rl2 + rl3)
w1 <- rl1/(rl0 + rl1 + rl2 + rl3)
w2 <- rl2/(rl0 + rl1 + rl2 + rl3)
w3 <- rl3/(rl0 + rl1 + rl2 + rl3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make some plots to visualize
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(y ~ s, data = data, ylab = 'Count (y)', xlab = 'Soil richness',
     col = b + 1, cex.lab = 2, pch = 19)
points(exp(p0) ~ s, data = data, cex = 3)


plot(y ~ s, data = data, ylab = 'Count (y)', xlab = 'Soil richness',
     col = b + 1, cex.lab = 2, pch = 19)
points(exp(p1) ~ s, data = data, cex = 3)

plot(y ~ jitter(b), data = data, ylab = 'Count (y)', xlab = 'Soil richness',
     col = b + 1, cex.lab = 2, pch = 19, xaxt = 'n')
axis(side = 1, at = c(0,1), labels = c('Unburned', 'Burned'))
points(exp(p2) ~ b, data = data, cex = 3, col = b+1)

plot(y ~ s, data = data, ylab = 'Count (y)', xlab = 'Soil richness',
     col = b + 1, cex.lab = 2, pch = 19)
points(exp(p3) ~ s, data = data, cex = 3, col = b+1)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's try a 4th hypothesis, perhaps there is an interaction b/w soil and fire
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(m4 <- glm(y ~ s*b, data = data, family = 'poisson'))
model.sel(m0,m1,m2,m3,m4)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note that m4 is almost exactly 2 AICc worse the model 3!
# it's simply m3, plus a junk parameter.
# It is NOT competitive with m3.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~











# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Change one of alpha[2] or alpha[3] parameters (used to set effects of soil (alpha[2]))
#    and fire (alpha[3]) to 0... how does that affect model selection?
#    Why?! :)
# 2) which of those four models are competitive? Should we remove one from the model set?
#    Think carefully about the parameter penalty and what it means.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's look at a simpler example
# log-likelihoood for normal distribution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmpn <- 100
tmpx <- rnorm(tmpn, 0, 1)
tmpy <- rnorm(tmpn, tmpx, 1)
plot(tmpy ~ tmpx, ylab = 'Response', xlab = 'Covariate', cex.lab = 2)
summary(tmpm0 <- glm(tmpy ~ 1))
summary(tmpm1 <- glm(tmpy ~ tmpx))

plot(tmpy ~ tmpx, ylab = 'Response', xlab = 'Covariate', cex.lab = 2)
points(tmpm0$fitted.values ~ tmpx, cex = 3)

plot(tmpy ~ tmpx, ylab = 'Response', xlab = 'Covariate', cex.lab = 2)
points(tmpm1$fitted.values ~ tmpx, cex = 3)


summary(tmpm0)
# calculate log-likelihood and AIC for null model
tmpl0 <- -(tmpn/2)*log(2*pi) - (tmpn/2)*log(sigma(tmpm0)^2) -
  (1/(2 * sigma(tmpm0)^2)) * sum((tmpy - tmpm0$fitted.values)^2)
2*2 + -2 * tmpl0

# calculate log-likelihood and AIC for true model
tmpl1 <- -(tmpn/2)*log(2*pi) - (tmpn/2)*log(sigma(tmpm1)^2) -
  (1/(2 * sigma(tmpm1)^2)) * sum((tmpy - tmpm1$fitted.values)^2)
2*3 + -2 * tmpl1

rm(tmpn, tmpx, tmpy, tmpm0, tmpm1, tmpl0, tmpl1)

