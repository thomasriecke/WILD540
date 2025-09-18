# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(terra)
# library(geoR)
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
# points(c(25,25,49) ~ c(1,25,49))
head(g)






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a new box... here we'll use a function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nubox <- function(ncells, x, y){
  
  # error trap and message to make sure number of cells is even 
  # so it can be symmetrical around centroid
  if (ncells %% 2 != 0){
    stop("ncells is not an even number")
    geterrmessage()
  }
  
  # error trap and message to make sure the box fits on the landscape
  if ((x + ncells/2) > ng | (y + ncells/2) > ng | 
      (x - ncells/2) < 0 | (y - ncells/2) < 0){
    stop("You can't put that box there")
    geterrmessage()
  }
  
 box <- st_polygon(list(cbind(c(x-(ncells/2),x+(ncells/2),x+(ncells/2),x-(ncells/2),x-(ncells/2)),
                              c(y-(ncells/2),y-(ncells/2),y+(ncells/2),y+(ncells/2),y-(ncells/2)))))
 
 return(box = box)
}  
treat <- nubox(10, 30 ,30)
plot(treat, col = 'red', add = T)




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


set.seed(29)
epsilon <- rmvnorm(1, mu = 0, Sigma)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# random noise across landscape
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- st_as_sf(g)
tmp <- st_contains(treat,g)
tmp <- as.numeric(tmp[[1]])                     # these are the cells inside the treated area
z$v <- as.numeric(epsilon)


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
# add a treatment effect
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta <- 2
z$v[tmp] <- epsilon[tmp] + beta

plot(z, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# which cells are not in the treated area?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`%!in%` = Negate(`%in%`)
out <- which(c(1:nc) %!in% tmp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll sample the treated area (st) and the rest of the study area (su)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n1 <- 20
n2 <- 20
s1 <- sample(tmp, n1, replace = F)
s2 <- sample(out, n2, replace = F)
y1 <- z$v[s1]
y2 <- z$v[s2]

plot(z, 
     main = NULL,                           # title for plot (NULL to remove)
     breaks = seq(-5,5,length.out = 40),    # a numeric vector with the actual breaks
     key.pos = 4,                           # numeric; side to plot a color key: 1 bottom, 2 left, 3 top, 4 right
     reset = F)                             # if FALSE, keep the plot in a mode that allows adding further map elements 
points(p[s1], pch = 19, col = 'black')
points(p[s2], pch = 19, col = 'dodgerblue3')

par(mfrow = c(1,1), mar = c(5,5,2,2))
boxplot(y1, y2, names = c('Treated', 'Untreated'), col = c('grey10','dodgerblue3'), boxwex = 0.5)

# ?t.test()
t.test(y1,y2)
t = (mean(y1) - mean(y2))/sqrt((var(y1)/n1) + (var(y2)/n2))
# if t > 2, p < 0.05
sig <- t > 2


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 1: Take a minute and think about some potential issues with 
#            this approach, i.e., simply comparing an area we did
#            something to with an area we didn't do anything to
# 
#            maybe one of our simulations highlighted this issue already :)
#            How could we address this? We'll do this soon in class!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Problem 2:
# proving that t-tests work. here we simulate two datasets using
# identical parameters (mean 0, st. dev. 1) and then compare them
# how often do you expect them to 'significantly' differ by random chance?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 20
nSims <- 10000
p <- NULL
for (i in 1:nSims){
  y <- rnorm(n, 0, 1)
  x <- rnorm(n, 0, 1)
  p[i] <- t.test(y,x)$p.value
}
# 95% of the time it says there is no difference
# x and y were both randomly drawn from normal(0,1)
length(which(p > 0.05))/nSims


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 3:
# change the effect size (beta)? what changes
# now change the sample size (n)?
# think about adding a second loop to try different sample sizes with
# different sample sizes...?
# different effect sizes...?
# different significance levels? (can extract post-hoc...)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 20
nSims <- 10000
p <- NULL
beta <- 1             # effect size
alpha <- 0.05         # alpha value or 'significance level'
for (i in 1:nSims){
  y <- rnorm(n, beta, 1)
  x <- rnorm(n, 0, 1)
  p[i] <- t.test(y,x)$p.value
}

length(which(p > 0.05))/nSims  # statistical power





