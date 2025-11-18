# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script includes an example of SEM/occupancy models:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(latex2exp)
library(reshape2)
library(jagsUI)
library(nimble)

set.seed(1234)
n <- 300
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1a) Let's simulate proportion aspen across n sites as a z-standardized covariate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aspen <- rnorm(n, 0, 1)
hist(aspen, breaks = 50, main = NULL, xlab = 'Z-standardized proportion aspen (x)', cex.lab =2)
# large positive values indicate lots of aspen and vice versa 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1b) Let's simulate latent occupancy as a function of proportion aspen
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(0,2)
psi <- plogis(alpha[1] + alpha[2] * aspen)
z <- rbinom(n, 1, psi)

par(mfrow = c(2,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
plot(psi ~ aspen, ylab = TeX("Probability of occupancy ($\\psi$)"),
     xlab = 'Z-standardized proportion aspen (x)', las = 1, cex.lab = 1.5)
plot(jitter(z,0.1) ~ aspen, ylab = TeX("Site occupancy (z)"),
     xlab = 'Z-standardized proportion aspen (x)', las = 1, cex.lab = 1.5)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1b) Let's simulate our observations as a function of latent occupancy (z)
#     and occasion-specific detection probability
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
J <- 4
beta <- c(-1,0,0,-1)
p <- plogis(beta)
y <- matrix(0, n, J)

for (j in 1:J){
  y[,j] <- rbinom(n, z, p[j])
}

par(mfrow = c(2,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
plot(p ~ c(1:J), ylab = TeX("Detection probability ($\\p$)"),
     xlab = 'Occasion', las = 1, cex.lab = 1.5, type = 'b', ylim = c(0,1))
barplot(colSums(y) ~ c(1:4), ylab = TeX("Detections (y)"),
     xlab = 'Occasion', las = 1, cex.lab = 1.5)


par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
plot(p ~ c(1:J), ylab = TeX("Detection probability ($\\p$)"),
     xlab = 'Occasion', las = 1, cex.lab = 1.5, type = 'b', ylim = c(0,1))


1 - prod(1-p)
p


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Let's format our data for Nimble and write a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nimble model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  alpha[1] ~ dlogis(0, 1)
  alpha[2] ~ dnorm(0,1)
  
  for (j in 1:J){
    beta[j] ~ dlogis(0,1)
    logit(p[j]) <- beta[j]
  }
  
  for (i in 1:n){
    
    logit(psi[i]) <- alpha[1] + alpha[2] * x[i]
    z[i] ~ dbern(psi[i])
    
    for (j in 1:J){
      y[i,j] ~ dbern(z[i] * p[j])
    }

  }
  
} )




# Bundle data (same as for JAGS)
z <- rowSums(y)
z[z > 0] <- 1
z[z == 0] <- NA
str(dataList <- list(y = y, x = aspen, n = n, z = z, J = J)) # not shown




init.z <- z
init.z[is.na(z)] <- 0
# Can use same function to generate starting values as for JAGS
inits <- function()  list (alpha = c(0,2), beta = rep(0,J), z = init.z)

# Parameters monitored: same as before
params <- c("alpha","beta")

# MCMC settings
ni <- 10000 ; nb <- 5000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
  m <- nimbleMCMC(code = model, constants = dataList,
                  inits = inits, monitors = params, nburnin = nb, niter = ni,
                  thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(m, head) # print first 5 values of each chain (not shown)
samps <- do.call(rbind.data.frame, m$samples)
str(samps)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
#
# 1a) Change the detection probabilities (work with a partner to increase or decrease
#     them). How does that affect your estimates of occupancy parametesr? Why?
# 
# 1b) Change the number of visits to 6. Note that you'll have to add a couple of betas
#     so that there is a detection probability for every occasion...
#     What happens to your estimates of alpha (occupancy parameters), 
#     specifically the uncertainty associated with those estiamtes.  Why?
#
# 2) work through plotting code of expected ruffed grouse occupancy 
#    regressed against aspen. Make sure it makes sense. Note how the plot changes
#    as you change detection probabilities
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

res <- 100
px <- seq(-3,3,length.out = res)
py <- matrix(NA, nrow(samps), res)
qy <- matrix(NA, 5, res)
for (j in 1:res){
  py[,j] <- plogis(samps[,1] + samps[,2] * px[j])
  qy[,j] <- quantile(py[,j], c(0.025,0.25,0.5,0.75,0.975))
}

plot(qy[3,] ~ px, type = 'l', lwd = 3,
     ylim = c(0,1),
     ylab = TeX("Probability of occupancy ($\\psi$)"),
     cex.lab = 2,
     xlab = 'Z-standardized proportion aspen (x)')
lines(qy[1,] ~ px, lty = 2, lwd = 2)
lines(qy[5,] ~ px, lty = 2, lwd = 2)
lines(qy[2,] ~ px, lty = 5, lwd = 1)
lines(qy[4,] ~ px, lty = 5, lwd = 1)






