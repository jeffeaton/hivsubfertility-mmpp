################################################################################
####                                                                        ####
####  This script fits the MMPP model.                                      ####
####                                                                        ####
####  Created by Jeff Eaton on 17 January 2016.                             ####
####                                                                        ####
################################################################################


#####################
####  Fit model  ####
#####################

setwd("~/Documents/Research/hiv-subfertility/mmpp/")
load("alpha-mmpp-data.RData")

source("mmpp-model.R")

par2 <- c(-4.81089224, -3.99204589, -3.89437682, -4.08546074, -4.36714715, -4.21024969, -7.29942232, -2.04898019, -1.32786472, -1.38946660, -1.51588339, -1.82244064, -2.45752807, -3.85530118, 1.00242190, 0.08917528, -0.04615499, -0.07043248, -0.06126176, -0.21694795, 0.20863540, -0.29060551, -2.01196198, -2.70513295)



system.time(print(ll(par2, epis.hivp, dat.hivn, epis.entry)))

options(width=500)

fit.alpha <- optim(par2, ll, epis.hivp = epis.hivp, dat.hivn = dat.hivn, epis.entry = epis.entry,
                   method="BFGS", control=list(fnscale=-1, trace=4, REPORT=1), hessian=TRUE)

fit.alpha2 <- optim(rep(log(0.3), 24), ll, epis.hivp = epis.hivp, dat.hivn = dat.hivn, epis.entry = epis.entry,
                   method="BFGS", control=list(fnscale=-1, trace=4, REPORT=1), hessian=TRUE)

save(fit.alpha, fit.alpha2, file="alpha-mmpp-optim-fit_2015-01-17.RData")

nll <- function(...) -ll(...)

nlm.alpha <- nlm(nll, par2, epis.hivp = epis.hivp, dat.hivn = dat.hivn, epis.entry = epis.entry,
                 fscale=-1, hessian=TRUE, print.level=2)

save(nlm.alpha, nlm.alpha2, file="alpha-mmpp-nlm-fit_2015-01-17.RData")

nlm.alpha2 <- nlm(nll, rep(log(0.3), 24), epis.hivp = epis.hivp, dat.hivn = dat.hivn, epis.entry = epis.entry,
                 fscale=-1, hessian=TRUE, print.level=2)
