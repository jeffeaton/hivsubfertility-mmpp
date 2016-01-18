################################################################################
####                                                                        ####
####  Likelihood for MMPP model for fertility by stage of HIV infection,    ####
####  with left truncation.                                                 ####
####                                                                        ####
####  Created by Jeff Eaton on 4 January 2016.                              ####
####                                                                        ####
################################################################################


## Fixed parameters
mu.hivn <- 0.005
mu <- c(0, 0.004, 0.010, 0.026, 0.435) + mu.hivn
pi0 <- c(0.58, 0.23, 0.16, 0.03)
sigma <- 1/c(6.37, 2.86, 3.54)
NS <- length(mu)


eAt.eigdec <- function(t, eigval, V, Vinv){
  V %*% (exp(t*eigval) * Vinv)
}

createQobj <- function(alpha, rho, lambda, sigma, pi0, mu){

  ## alpha: age specific fertility rate
  ## rho: stage-specific subfertility
  ## lambda: HIV incidence rate
  ## sigma: stage progression rates
  ## pi0: initial stage distribution after seroconversion
  ## mu: 

  NS <- length(mu)  # uninfected + disease stages
  
  Q  <- diag(-c(lambda, sigma, 0) - alpha*c(1, rho) - mu)
  Q[1, 2:NS] <- lambda*pi0
  Q[NS*(2:(NS-1)) + 2:(NS-1)] <- sigma  # upper off-diagnonal from 3rd column

  Qe <- eigen(Q, symmetric=FALSE)

  Emat <- list(birth    = diag(alpha * c(1, rho)),
               death    = diag(mu),
               lastneg  = diag(c(1, rep(0, NS-1))),
               firstpos = diag(c(0, rep(1, NS-1))),
               cens     = diag(NS))

  list(Q       = Q,
       eigvalQ = Qe$values,
       V       = Qe$vectors,
       Vinv    = solve(Qe$vectors),
       Emat    = Emat)
}


ll <- function(par, epis.hivp, dat.hivn, epis.entry){

  print(round(par, 8))

  lambda <- exp(par[1:7]) # incidence rate for each age group
  alpha <- exp(par[8:14]) # age-specfic fertility
  frr.age <- exp(par[15:21])
  frr.stage <- c(1, exp(par[22:24]))  # subfertility by stage
  
  if(any(exp(par) == Inf) | any(alpha %o% frr.stage == Inf))
    return(-Inf)

  Qlist <- lapply(1:nrow(Qlist.mm),
                  function(i){
                    aidx <- Qlist.mm[i,"agegr"]
                    createQobj(ifelse(Qlist.mm[i,"fertelig"], alpha[aidx], 0),
                               frr.age[aidx]*frr.stage, lambda[aidx], sigma, pi0, mu)})
  
  ll.hivn(dat.hivn, Qlist) + sum(sapply(epis.hivp, ll.hivp, Qlist)) - log.pentry(epis.entry, Qlist)
}

ll.hivn <- function(dat.hivn, Qlist){
  q.hivn <- sapply(Qlist, function(x) x$Q[1])
  logrho.hivn <- log(sapply(Qlist, function(x) x$Emat$birth[1]))
  logmu.hivn <- log(sapply(Qlist, function(x) x$Emat$death[1]))
  ## logmu.hivn <- replace(logmu.hivn, logmu.hivn == -Inf, -1e6)

  ll.hivn <- sum(q.hivn[dat.hivn$Qidx] * dat.hivn$dur) +
    sum(replace(logrho.hivn[dat.hivn$Qidx] * dat.hivn$births, dat.hivn$births==0, 0)) +
    sum(logmu.hivn[dat.hivn$Qidx] * dat.hivn$deaths)

  return(ll.hivn)
}

ll.hivp <- function(e, Qlist){

  Qidx <- e[["Qidx"]]
  dur <- e[["dur"]]
  type <- e[["type"]]

  forwp.t <- c(1, rep(0, NS-1))  # forward probabilities at time t.
  for(ii in 1:nrow(e)){
    Qobj <- Qlist[[Qidx[ii]]]
    forwp.t <- forwp.t %*%
      eAt.eigdec(dur[ii], Qobj[["eigvalQ"]], Qobj[["V"]], Qobj[["Vinv"]]) %*%
      Qobj[["Emat"]][[type[ii]]]
  }
  return(log(sum(forwp.t)))
}

log.pentry <- function(epis.entry, Qlist){

  Qidx <- epis.entry[["Qidx"]]
  dur <- epis.entry[["dur"]]
  entry <- epis.entry[["entry"]]
  
  psurv <- c(1, rep(0, NS-1))
  lpentry <- 0
  for(ii in 1:nrow(epis.entry)){
    Qobj <- Qlist[[Qidx[ii]]]
    psurv <- psurv %*% eAt.eigdec(dur[ii], Qobj[["eigvalQ"]], Qobj[["V"]], Qobj[["Vinv"]])
    if(entry[ii])
      lpentry <- lpentry + log(sum(psurv))
  }
  return(lpentry)
}
