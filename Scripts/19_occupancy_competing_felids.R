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
# simulate our two covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
forest <- rnorm(n, 0, 1)
snow <- rnorm(n, 0, 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate bobcat occupancy probability (psiB)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha0 <- -1
alpha1 <- -3
psiB <- plogis(alpha0 + alpha1 * snow)
plot(psiB ~ snow) # bobcats aren't present in deep snow

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate bobcat occupancy (zB)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zB <- rbinom(n, 1, psiB)
table(zB)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate lynx occupancy probability (psiL)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta0 <- 0
beta1 <- 0.5
beta2 <- -2
psiL <- plogis(beta0 + beta1 * forest + beta2 * zB)
plot(psiL ~ forest) 
# these two lines of points represent lynx occupancy probabilities at sites
# with (upper) and without (lower) bocats

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate lynx occupancy (zL)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zL <- rbinom(n, 1, psiL)
table(zL)

table(zL,zB)
# 97 sites have neither species, only 19 have both


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate detection data for bobcat (yB) and lynx (yL)
# we'll sample each site for four weeks (K)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
K <- 4
yB <- rep(0, n)
yL <- rep(0, n)
p <- 0.5
for (i in 1:n){
  yB[i] <- rbinom(1, K * zB[i], p)
  yL[i] <- rbinom(1, K * zL[i], p)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(yB,zB)
table(yL,zL) 





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Let's format our data for Nimble and write a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nimble model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- nimbleCode( {
  
  alpha0 ~ dlogis(0, 1)
  alpha1 ~ dnorm(0, 0.1)
  
  beta0 ~ dlogis(0, 1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  pB ~ dbeta(1, 1)
  pL ~ dbeta(1, 1)
  
  for (i in 1:n){
    
    logit(psiB[i]) <- alpha0 + alpha1 * s[i]
    logit(psiL[i]) <- beta0 + beta1 * f[i] + beta2 * zB[i]
    # note we use zB as a covariate instead of psiB because we care 
    # about the effect of bobcat presence (zB)...
    # not the effect of bocat probability of presence (psiB)
    
    zB[i] ~ dbern(psiB[i])
    zL[i] ~ dbern(psiL[i])
    
    yB[i] ~ dbin(pB, zB[i] * K)
    yL[i] ~ dbin(pL, zL[i] * K)
    
  }
  
} )




# Bundle data (same as for JAGS)
str(dataList <- list(yB = yB, yL = yL, 
                     f = forest,
                     s = snow,
                     n = n,
                     K = K)) # not shown





# Can use same function to generate starting values as for JAGS
inits <- function() list ()

# Parameters monitored: same as before
params <- c("alpha0","alpha1",
            "beta0","beta1","beta2",
            'pB','pL')

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


boxplot(samps)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK (don't turn this in, just think about it)
#
# 1) Make a plot of the effect of bobcat occupancy on lynx occupancy
#
# 2) Work with a neighbor to think about a different model parameterization given
#    alternative hypotheses
#    
# 3) revise the code above the model to strengthen or weaken the effects
#    of bobcat on lynx. What happens to the raw 'contingency tables'
#    i.e., table(zB,zL)
#    Why?
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n.iter <- nrow(samps)
res <- 100
xf <- seq(-3,3,length.out = res)
xs <- seq(-3,3,length.out = res)
psiB <- matrix(NA, n.iter, res)
psiL1 <- matrix(NA, n.iter, res)
psiL0 <- matrix(NA, n.iter, res)
qpsiB <- matrix(NA, res, 5)
qpsiL1 <- matrix(NA, res, 5)
qpsiL0 <- matrix(NA, res, 5)

for (j in 1:res){
  psiB[,j] <- plogis(samps$alpha0 + samps$alpha1 * xs[j])
  psiL1[,j] <- plogis(samps$beta0 + samps$beta1 * xf[j] + samps$beta2)  
  psiL0[,j] <- plogis(samps$beta0 + samps$beta1 * xf[j])
  qpsiB[j,] <- quantile(psiB[,j], c(0.025,0.05,0.5,0.95,0.975))
  qpsiL1[j,] <- quantile(psiL1[,j], c(0.025,0.05,0.5,0.95,0.975)) 
  qpsiL0[j,] <- quantile(psiL0[,j], c(0.025,0.05,0.5,0.95,0.975))   
}

mpsiL1 <- melt(psiL1)
mpsiL0 <- melt(psiL0)
mpsiB <- melt(psiB)

par(mfrow = c(1,2))
smoothScatter(mpsiL1$value ~ xf[mpsiL1$Var2], nrpoints = 0, ylim = c(0,1),
              xlab = 'Forest', ylab = 'P(Lynx occupancy | bobcat presence)')
lines(qpsiL1[,3] ~ xf, col = 'white', lwd = 3)
lines(qpsiL1[,1] ~ xf, col = 'white', lty = 2, lwd = 3)
lines(qpsiL1[,5] ~ xf, col = 'white', lty = 2, lwd = 3)

smoothScatter(mpsiL0$value ~ xf[mpsiL1$Var2], nrpoints = 0, ylim = c(0,1),
              xlab = 'Forest', ylab = 'P(Lynx occupancy | bobcat absence)')
lines(qpsiL0[,3] ~ xf, col = 'white', lwd = 3)
lines(qpsiL0[,1] ~ xf, col = 'white', lty = 2, lwd = 3)
lines(qpsiL0[,5] ~ xf, col = 'white', lty = 2, lwd = 3)



withoutB <- plogis(samps$beta0)
withB <- plogis(samps$beta0 + samps$beta2)
qwith <- quantile(withB, c(0.025,0.5,0.975))
qwithout <- quantile(withoutB, c(0.025,0.5,0.975))
par(mfrow = c(1,1))
vioplot::vioplot(withoutB, withB, drawRect = F,
        names = c('no bobcats','bobcats'))

arrows(1, qwithout[1], 1, qwithout[3], length = 0, col = 'white', lwd = 2)
arrows(2, qwith[1], 2, qwith[3], length = 0, col = 'white', lwd = 2)
points(c(qwithout[2],qwith[2]) ~ c(1:2), pch = 21, bg = 'white', cex = 3)

delta <- withB - withoutB
hist(delta, breaks = 100)

