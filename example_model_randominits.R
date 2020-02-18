#=================================================================
# File:
#   Example code for estimating a latent class joint model
#   using "random" initial values (see the paper's
#   Supplementary Materials for more details on the approach).
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

#----------

# Function to draw initial values
getinits <- function() {
  c(
    intercept_class1 = rnorm(1, 0, 1),
    intercept_class2 = rnorm(1, 0, 1),
    intercept_class3 = rnorm(1, 0, 1),
    intercept_class4 = rnorm(1, 0, 1),
    event1_weibull1_class1 = rnormsq(1, 0, 1),
    event1_weibull2_class1 = rnormsq(1, 0, 1),
    event1_weibull1_class2 = rnormsq(1, 0, 1),
    event1_weibull2_class2 = rnormsq(1, 0, 1),
    event1_weibull1_class3 = rnormsq(1, 0, 1),
    event1_weibull2_class3 = rnormsq(1, 0, 1),
    event1_weibull1_class4 = rnormsq(1, 0, 1),
    event1_weibull2_class4 = rnormsq(1, 0, 1),
    event1_weibull1_class5 = rnormsq(1, 0, 1),
    event1_weibull2_class5 = rnormsq(1, 0, 1),
    event2_weibull1_class1 = rnormsq(1, 0, 1),
    event2_weibull2_class1 = rnormsq(1, 0, 1),
    event2_weibull1_class2 = rnormsq(1, 0, 1),
    event2_weibull2_class2 = rnormsq(1, 0, 1),
    event2_weibull1_class3 = rnormsq(1, 0, 1),
    event2_weibull2_class3 = rnormsq(1, 0, 1),
    event2_weibull1_class4 = rnormsq(1, 0, 1),
    event2_weibull2_class4 = rnormsq(1, 0, 1),
    event2_weibull1_class5 = rnormsq(1, 0, 1),
    event2_weibull2_class5 = rnormsq(1, 0, 1),
    rxfirstdatecat10to14_event1 = rnorm(1, 0, 0.25),
    rxfirstdatecat10to14_event2 = rnorm(1, 0, 0.25),
    femaleYes_event1            = rnorm(1, 0, 0.25),
    femaleYes_event2            = rnorm(1, 0, 0.25),
    lungynYes_event1            = rnorm(1, 0, 0.25),
    lungynYes_event2            = rnorm(1, 0, 0.25),
    coronaryynrevNo_event1      = rnorm(1, 0, 0.25),
    coronaryynrevNo_event2      = rnorm(1, 0, 0.25),
    cvdynrevNo_event1           = rnorm(1, 0, 0.25),
    cvdynrevNo_event2           = rnorm(1, 0, 0.25),
    pvdynYes_event1             = rnorm(1, 0, 0.25),
    pvdynYes_event2             = rnorm(1, 0, 0.25),
    diabetesynrevNo_event1      = rnorm(1, 0, 0.25),
    diabetesynrevNo_event2      = rnorm(1, 0, 0.25),
    prdcatGN_event1             = rnorm(1, 0, 0.25),
    prdcatGN_event2             = rnorm(1, 0, 0.25),
    prdcatHypertension_event1   = rnorm(1, 0, 0.25),
    prdcatHypertension_event2   = rnorm(1, 0, 0.25),
    prdcatOtherUncertain_event1 = rnorm(1, 0, 0.25),
    prdcatOtherUncertain_event2 = rnorm(1, 0, 0.25),
    racecat2Aboriginal_event1   = rnorm(1, 0, 0.25),
    racecat2Aboriginal_event2   = rnorm(1, 0, 0.25),
    racecat2Asian_event1        = rnorm(1, 0, 0.25),
    racecat2Asian_event2        = rnorm(1, 0, 0.25),
    racecat2Maori_event1        = rnorm(1, 0, 0.25),
    racecat2Maori_event2        = rnorm(1, 0, 0.25),
    agerrt50_event1             = rnorm(1, 0, 0.25),
    agerrt50_event2             = rnorm(1, 0, 0.25),
    intercept_class1 = runif(1, 0, 1),
    intercept_class2 = runif(1, 0, 1),
    intercept_class3 = runif(1, 0, 1),
    intercept_class4 = runif(1, 0, 1),
    intercept_class5 = runif(1, 0, 1),
    ns1_class1 = rnorm(1, 0, 0.25),
    ns1_class2 = rnorm(1, 0, 0.25),
    ns1_class3 = rnorm(1, 0, 0.25),
    ns1_class4 = rnorm(1, 0, 0.25),
    ns1_class5 = rnorm(1, 0, 0.25),
    ns2_class1 = rnorm(1, 0, 0.25),
    ns2_class2 = rnorm(1, 0, 0.25),
    ns2_class3 = rnorm(1, 0, 0.25),
    ns2_class4 = rnorm(1, 0, 0.25),
    ns2_class5 = rnorm(1, 0, 0.25),
    ns3_class1 = rnorm(1, 0, 0.25),
    ns3_class2 = rnorm(1, 0, 0.25),
    ns3_class3 = rnorm(1, 0, 0.25),
    ns3_class4 = rnorm(1, 0, 0.25),
    ns3_class5 = rnorm(1, 0, 0.25),
    varcov1 = runif(1, 0, 1),
    stderr  = runif(1, 0, 1)
  )
}

#----------

# Fit a latent class joint model with 5 classes, using
# a randomly drawn set of initial values
seed <- rint()       # generate a random seed
set.seed(seed)       # set the seed for reproducibility
inits <- getinits()  # randomly draw a set of initial values
mod <- Jointlcmm(    # fit the latent class joint model
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
  B = inits,
  maxiter = 200)
mod <- structure(mod, inits = inits, seed = seed)
saveRDS(mod, "C:/replace/file/pathway/here.rds")

