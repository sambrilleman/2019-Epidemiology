#=================================================================
# File:
#   Example code for estimating a latent class joint model
#   using our "Strategy B" for specifying initial values (see the
#   paper's Supplementary Materials for more details on the approach).
#
# Paper:
#   "Longitudinal changes in body mass index and the competing
#   outcomes of death and transplant in patients undergoing
#   hemodialysis: a latent class joint modelling approach"
#   by Brilleman et al.
#
# Copyright: Sam Brilleman, 2018
#=================================================================

library(lcmm)
library(splines)

#----------

# Function to draw a random integer (used as the seed)
rint <- function() {
  as.integer(runif(1, 1, .Machine$integer.max))
}

# Set seed
seed <- rint() # generate a random seed
set.seed(seed) # set the seed for reproducibility

#----------

# Fit a latent class joint model with only 1 class
# that can be used to obtain initial values for the
# latent class joint model with >1 classes
initmod <- Jointlcmm(
  fixed = bmi ~ ns(t0, 3),
  random = ~ 1,
  subject = "randid",
  ng = 1,
  data = data,
  survival = Surv(competT, competD) ~
    cause(rxfirstdatecat) +
    cause(female) +
    cause(lungyn) +
    cause(coronaryynrev) +
    cause(cvdynrev) +
    cause(pvdyn) +
    cause(diabetesynrev) +
    cause(prdcat) +
    cause(racecat2) +
    cause(agerrt50),
  logscale = FALSE,
  maxiter = 200)

#----------

# Fit a latent class joint model with 5 classes, using
# the 1 class model to generate initial values
mod <- gridsearch(
  Jointlcmm(
    fixed = bmi2 ~ ns(t0, 3),
    mixture = ~ ns(t0, 3),
    random = ~ 1,
    subject = "randid",
    ng = 5,
    data = data,
    survival = Surv(competT, competD) ~
      cause(rxfirstdatecat) +
      cause(female) +
      cause(lungyn) +
      cause(coronaryynrev) +
      cause(cvdynrev) +
      cause(pvdyn) +
      cause(diabetesynrev) +
      cause(prdcat) +
      cause(racecat2) +
      cause(agerrt50),
    logscale = FALSE,
    B = initmod,
    maxiter = 200),
  rep = 3, maxiter = 50, minit = initmod)

#---------

# If the model failed to converge, then rerun using the
# initial values for the longitudinal submodel that are
# obtained from a latent class mixed model with five
# classes but no survival component.

# Function to rerun the joint model using initial values
# from a separate longitudinal model fit using lcmm::hlme
rerun_using_hlmeinits <- function(mod) {
  mc1 <- mc2 <- mod$call
  mc1[[1]] <- hlme
  mc1$survival <- NULL
  mc1$logscale <- NULL
  mc1$B <- NULL
  hlmemod <- eval(mc1)
  hlmeinits <- hlmemod$best
  a <- 1
  b <- hlmemod$ng - 1
  c <- hlmemod$ng
  d <- length(hlmeinits)
  hlmeinits1 <- hlmeinits[a:b]
  hlmeinits2 <- hlmeinits[c:d]
  inits <- mod$best
  a <- 1
  b <- mod$ng - 1
  c <- length(inits) - length(hlmeinits2) + 1
  d <- length(inits)
  inits[a:b] <- hlmeinits1
  inits[c:d] <- hlmeinits2
  mc2$B <- inits
  mod <- eval(mc2)
  return(mod)
}

# Then rerun the model if required
if (mod$conv == 1) {
  mod <- structure(mod, rerun = FALSE, seed = seed)
} else {
  mod <- rerun_using_hlmeinits(mod)
  mod <- structure(mod, rerun = TRUE, seed = seed)
}
saveRDS(mod, "C:/replace/file/pathway/here.rds")

