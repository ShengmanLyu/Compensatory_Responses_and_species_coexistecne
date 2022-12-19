#*****************************************************************************
# Functions ----
# Author: Shengman Lyu
# E-mail:shengman.lyu@gmail.com
# Date updated: 19.12.2022
# More information can be found in Lyu, S. and Alexander, J. (2022) Compensatory responses of vital rates attenuate impacts of competition on population growth and promote coexistence
#*****************************************************************************
library(MuMIn)
library(car)
library(ggplot2)

# The codes provided here define functions used for 
#  - builing intergral projection models (IPM)
#  - basic analyses of IPMs: projection, population growth rates, LTRE
#  - coexistence analyses: coexistence outcomes, niche differences, relative fitness differences
#  - statistical analyses: different types of confidence intervals

#************************************************************
# ****** Functions for making intergral projection model (IPMs) ****** ----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vital rate parameters----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par <- c("surv.int", "surv.slope", 
         "growth.int", "growth.slope", "growth.sd", 
         "flowering.int", "flowering.slope", 
         "fecundity.int", "fecundity.slope",
         "germination.prob",
         "establishment.funiche.int", 
         "establishment.comp.int",
         "seedling.size.mean",
         "seedling.size.sd",
         "L", "U")
length(par)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vital rate functions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. probability of surviving (logit)
survival.z=function(z,params) {
  # linear predictor
  u <- params$surv.int+params$surv.slope*z
  # inverse logit transformation
  survival.p <- exp(u)/(1+exp(u))
  survival.p[u>700] <- 1 # exp(710) gives Inf values
  # return the probability of survival
  return(survival.p)
}

# 2. growth function (gaussian)
growth.z1z <- function(z1,z, params) {
  # mean size time t+1, z1
  mu <- params$growth.int+params$growth.slope*z
  # sd about mean
  sig <- params$growth.sd
  # probability density of new size z1 for current size = z
  pd <- dnorm(z1,mean=mu,sd=sig)
  # return probability density
  return(pd)
}

# 3 probaility of flowering (logit)
flowering.z <- function(z, params) {
  # linear predictor
  u <- params$flowering.int+params$flowering.slope*z
  # probability of flowering, inverse logit
  flowering.p <- exp(u)/(1+exp(u))
  flowering.p[u>700] <- 1 # exp(710) gives Inf values
  # return probability of flowering
  return(flowering.p) 
}

# 4. seed production (Gaussian and log)
seed.z.gaussian <- function(z,params) {
  # linear predictor
  u = params$fecundity.int+params$fecundity.slope*z
  seeds <- exp(u) 
  # nagetative seeds to 0
  seeds[seeds < 0] = 0 
  # return seed production
  return(seeds)
}

# 5 probability of seed germination (constant)
germination <- function(params) {
  # get the probability of seedling establishment
  germination.p <- params$germination.prob
  # return
  return(germination.p)
}

# 6.1 probability of seedling establishment in the absence of competition (logistic)
establishment.fun <- function(params) {
  # establishment using FuNiche data
  u.funiche <- params$establishment.funiche.int
  establishment.prob.funiche <- exp(u.funiche)/(1+exp(u.funiche))
  # return
  return(establishment.prob.funiche)
}

# 6.2 probability of seedling establishment in the presence of competition (logistic)
establishment.com <- function(params) {
  # establishment with background species
  u.comp <- params$establishment.comp.int
  establishment.prob.comp <- exp(u.comp)/(1+exp(u.comp))
  # return
  return(establishment.prob.comp)
}

# 7 seedling size (constant)
seedling.z1 <- function(z1, params) {
  # probability density of recruits
  rpd <- dnorm(z1, mean=params$seedling.size.mean, sd=params$seedling.size.sd)
  # return probability density
  return(rpd)
}

# G kernel: survival-growth
G.z1z <- function(z1,z,params) {
  # combine survival and growth
  sg <- survival.z(z, params) * growth.z1z(z1,z,params)
  # return
  return(sg)
}

# R kernel: reproduction
R.z1z=function(z1,z,params) {
  # calculate fecundity kernel
  f = flowering.z(z, params) * 
    seed.z.gaussian(z, params) * 
    germination(params) *
    establishment.fun(params) * establishment.com(params) *
    seedling.z1(z1, params)
  # return
  return(f)
}

# full kernel
full.z1z = function(z1,z,params) {
  # combine growth and reproduction kernel
  survival.z(z, params) * growth.z1z(z1,z,params) + 
    R.z1z(z1,z,params)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to make an IPM kernel ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To make an IPM kernel
mk.kernel <- function(L, U, par, n) {
# n: number of meshpoints
# par: vital rates
# L, U: lower and upper limit of size
# d: density of competitor
# fun: full kernel 
# mesh points 
  h <- (U - L)/n
  meshpts <- L + ((1:n) - 1/2) * h
  
  # kernel
  G <- h * (outer(meshpts, meshpts, G.z1z, params = par))
  R <- h * (outer(meshpts, meshpts, R.z1z, params = par))
  K <- G+R
  
  return(list(meshpts = meshpts, G=G, R=R, K = K))
}

#************************************************************
# ****** Functions for basic analyses for IPMs ******----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate population growth rate (lambda) by iteration----
# This is faster for 1000 x 1000 or bigger matrices
# Caution when the resulted lambda is big (> 50). maybe caused by iteractions!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda.iteration = function(k, tol=1e-8, i.max=10000) {
  # k: IPM full kerbel
  # tol: threshold for difference in iterated lambdas
  # i.max: maximum number of iterations
  qmax=10*tol
  i = 1 # number of interations
  lam=1
  x=rep(1,nrow(k))   
  k=Matrix(k)
  while(qmax>tol & i < i.max) {
    x1 = k%*%x
    qmax=sum(abs(x1-lam*x))  
    lam=sum(x1)
    i = i + 1
    # lambda is 0 when population is 0
    if(lam == 0) { x = 0; qmax = 0 }
    else x=x1/lam # one individual in total
  } 
  # lambda: converged population grwoth rate
  # w: stablized size distribution
  # n: number of interations
  return(list(lambda=lam,w=x/sum(x), n = i))
} 

# test
#k = ipm.i$K

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate population growth rate (lambda) using eigenvalue----
# This is faster than using iteration for 250 x 250 or smaller matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda.eigen <- function(k, only.lambda=TRUE) {
  # k: IPM full kerbel
  # only.lambda: whether or not calcaulte lamdab only or also stable size distribution
  eigen.k = eigen(k, only.values = only.lambda)
  lam = abs(eigen.k$values[1]) 
  if(only.lambda) return(list(lambda=lam, w=NA))
  else { 
    eigen.vec = eigen.k$vectors[,1]
    w = abs(eigen.vec)/sum(abs(eigen.vec))
    return(list(lambda=lam, w=w)) 
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to standardise relaitve contribution of vital rates ----
# This standardisation makes relative contribution compaable across species
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ltre.standardise <- function(x, absolute.value=TRUE) {
  # x: a vector of relative contribution, including positive and negtaive values, of vital rates to be standardised
  # absolute.value: standardization based on absolute value
  
  # Check if NA in x
  if(sum(is.na(x)) >0) warning("NAs in the relative contribution!")
  
  # absolute value
  if(absolute.value)  x.positive <- abs(x)
  else x.positive <- x
  
  # standardisation
  x.sum <- sum(x.positive)
  x.divided <- x.positive/x.sum
  x.standard <- ifelse(x< 0, (x.divided*-1), x.divided)
  
  return(x.standard)
}

# test
ltre.standardise(c(-1,1,2,-2))

#************************************************************
# ****** Functions for coexistence analyses ******----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate sensitivity of species using intrinsic and invasion growth rates ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sensitivity.igr <- function(pgr.intrinsic, pgr.invasion) {
  # pgr.intrinsic: intrinsic growth rate of the focal species
  # pgr.invasion: invasion growth rate of the focal species
  
  # log
  pgr.intrinsic.log <- log(pgr.intrinsic)
  pgr.invasion.log <- log(pgr.invasion)
  
  # calculate sensitivity
  if(pgr.intrinsic.log > 0) {
    sens <- 1- pgr.invasion.log/pgr.intrinsic.log
  }
  else if(pgr.intrinsic.log < 0) {
    sens <- pgr.invasion.log/pgr.intrinsic.log - 1
  }
  # output
  return(sens)
}

# when intrinic > 0
intrinsic = exp(0.5)
invasion = exp(0.5 + seq(-0.8, 0.8, length.out = 10))
sensitivity.igr(pgr.intrinsic = intrinsic, pgr.invasion = invasion)

# when intrinic < 0
intrinsic = exp(-0.5)
invasion = exp(-0.5 + seq(-0.8, 0.8, length.out = 10))
sensitivity.igr(pgr.intrinsic = intrinsic, pgr.invasion = invasion)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate sensitivity and coexistence outcome using invasion growth rate ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coex.igr <- function(sps1,sps2, igr10, igr20, igr12, igr21) {
  # igr10, igr20: intrinsic growth rate in the absence of competitors
  # igr11: invasion growth rate of species 1 invading species 2
  # igr12: invasion growth rates of species 2 invaing species 1
  
  # calculate sensitivity
  s1 <- sensitivity.igr(igr10, igr12); names(s1) <- sps1
  s2 <- sensitivity.igr(igr20, igr21); names(s2) <- sps2
  ss <- c(s1, s2)
  
  # log
  igr10 <- log(igr10)
  igr20 <- log(igr20)
  igr12 <- log(igr12)
  igr21 <- log(igr21)
  
  # coexsitence outcome independent of the intrinsic growth rates
  if(igr12 > 0 & igr21 >0) {
    outcome <- "Coexistence"
    winner <- "sps1_sps2"
  }
  else if(igr12 > 0 & igr21 < 0) { 
    outcome <- "Competitive exclusion"
    winner <- "sps1"
  }
  else if(igr12 < 0 & igr21 > 0) { 
    outcome <- "Competitive exclusion"
    winner <- "sps2"
  }
  else if(igr12 < 0 & igr21 < 0) { 
    outcome <- "Priority effect"
    winner <- "none"
  }
  
  # superior
  if(s1 > s2) superior <- "sps2"
  else superior <- "sps1"
  
  return(list(sensitivity = ss, outcome = outcome, winner = winner, superior = superior)) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate ND and FD and coexistence outcome ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coex.ndfd <- function(s12, s21) {
  # s12: sensitivity of species 1 invading species 2
  s <- c(s12, s21)
  
  # ND
  nd <- 1 - sqrt(s12*s21)
  names(nd) <- "ND"
  
  # FD
  # afd <- exp( sqrt( (log(s12/gm) + log(s21/gm))/2 ) )
  # afd <- exp(sqrt(mean(log(s)**2) - (mean(log(s)))**2))
  rfd <- sqrt(s21/s12)
  names(rfd)  <- "RFD"
  
  # superior
  if(rfd > 1) superior <- "sps1"
  else superior <- "sps2"
  
  # Coexistence outcomes
  # make sure FD always > 1
  if(rfd < 1) rfd2 <- 1/rfd
  else rfd2 <- rfd
  
  # coexistence metrics
  coex.metric <- 1/ (rfd2*(1-nd))
  names(coex.metric) <- "Coexistence metric"
  
  if(rfd2 <= 1/(1-nd)) { 
    outcome <- "Coexistence"
    winner <- "sp1_sps2"
    }
  else if(rfd > 1) {
    outcome <- "Competitive exclusion"
    winner <- "sps1"
    }
  else if(rfd < 1) {
    outcome <- "Competitive exclusion"
    winner <- "sps2"
  }
  return(list(ndfd = c(nd, rfd, coex.metric), outcome = outcome, winner = winner, superior = superior) )
}

# test
coex.ndfd(0.6,2)
coex.ndfd(2, 0.6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to plot ND, FD and coexistence outcome ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot ND and RFD (log-scale) frame indicating coexistence outcomes
coex.frame <- function(x1=-0.5,x2=0.9) {
  # x1, x2 lower and upper limits of ND, in which x1 < -0.001
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  p <- ggplot() + 
    geom_line(data=data.frame(x=x, y=y), aes(x=x, y=y), inherit.aes = FALSE) +
    geom_line(data=data.frame(x=x, y=-y), aes(x=x, y=y), inherit.aes = FALSE) +
    geom_line(data=data.frame(x=x, y=rep(0,length(x))), aes(x=x, y=y), inherit.aes = FALSE, linetype="dashed")
  return(p)
}
coex.frame()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add frame
frame.coex = function(x1 = -0.01, x2= 0.91) {
  # x1, x2 lower and upper limits of ND, in which x1 < -0.001
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  return(data.frame(x=x,y=y))
}
frame.coex()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add polygan
# polygan of coexistence area
polygan.coex = function(x1 = 0.0001, x2= 0.91) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(x2,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(log(1/(1-xx1)), log(1/(1-xx2)),-log(1/(1-xx3)))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.coex(), aes(x=x,y=y)) +
  geom_polygon()

polygan.coex.no = function(x1 = -1, x2=0) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = seq(x2,x1,length.out = 100)
  xx3 = rep(x1,100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(-xx1, xx2, seq(x1,-x1, length.out=100))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.coex.no(), aes(x=x,y=y)) +
  geom_polygon()

polygan.prio = function(x1 = -3, x2= -0.001) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(-0.001,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(log(1/(1-xx1)), log(1/(1-xx2)),-log(1/(1-xx3)))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.prio(), aes(x=x,y=y)) +
  geom_polygon()

polygan.prio.no = function(x1 = 0, x2=1) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(x2,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(xx1,seq(x2,-x2, length.out=100), -xx3)
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.prio.no(), aes(x=x,y=y)) +
  geom_polygon()

#************************************************************
# ****** Functions for statistical analyses ******----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compute confidence interval using linear model ----
# same as: Hmisc::mean_cl_normal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci <- function(x) {
  fit <- lm(x~1)
  c(mean(x), confint(fit)[1, ])
  }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compute SE-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_se <- function(x) {
  u = mean(x)
  se = sqrt(var(x)/length(x))
  data.frame(y=u, ymin=u-1.96*se, ymax=u+1.96*se)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compute SD-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_sd <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  data.frame(y=u, ymin=u-1.96*sd, ymax=u+1.96*sd)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compute mean, sd, cv, min, man and range ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_cv <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  cv = abs(sd/u)
  min = range(x, ...)[1]
  max = range(x, ...)[2]
  rg = max - min
  data.frame(n = length(x[!is.na(x)]), mean=u, sd=sd, cv=cv, range = rg, min=min, max=max)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compute mean, sd, cv, min, man and range ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_sd <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  min = u - sd
  max = u + sd
  data.frame(n = length(x[!is.na(x)]), mean=u, sd=sd, min=min, max=max)
}

mean_sd.ggplot <- function(x) {
  x <- stats::na.omit(x)
  u = mean(x)
  sd = sqrt(var(x))
  min = u - sd
  max = u + sd
  data_frame(y = u, ymin = min, ymax = max)
  #new_data_frame(list(y = u, ymin = min, ymax = max), n = 1)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate percentile-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_quantile <- function(x, ...) {
  u = mean(x, ...)
  x.min = quantile(x, probs = 0.025, ...)
  x.max = quantile(x, probs = 0.975, ...)
  data.frame(y=u, ymin=x.min, ymax=x.max) 
}

median_ci_quantile <- function(x, ...) {
  u = median(x, ...)
  x.min = quantile(x, probs = 0.025, ...)
  x.max = quantile(x, probs = 0.975, ...)
  data.frame(y=u, ymin=x.min, ymax=x.max) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate percentile-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
is.ci.significant <- function(x.min, x.max) {
  if((x.min > 0) == (x.max > 0)) is.sig <- TRUE
  else is.sig <- FALSE
  is.sig
}