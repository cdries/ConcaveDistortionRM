library(nleqslv)

############################
### Distortion functions ###
############################

# We use the convention of inserting zero whenever the input is outside the
# domain of the derivative of the distortion function [0, 1].
# In addition, we put the inverse to zero outside it's domain, as defined as
# the range of the derivative of the distortion function in [Ø, 1].

### Power distortion risk measure - p in (0, 1)
# distortion function
gPower <- function(x, p = 0.5) {
  y <- x^p
  y[x >= 1] <- 1
  y[x <= 0] <- 0
  return (y)
}

# derivative of distortion function
DgPower <- function(x, p = 0.5) {
  y <- p * x^(p - 1)
  y[x > 1] <- 0
  y[x <= 0] <- 0
  return (y)
}

# inverse of derivative of distortion function
DgPowerInv <- function(x, p = 0.5) {
  y <- rep(0, length(x))
  ind_dom <- which(x >= p)
  y[ind_dom] <- (x[ind_dom] / p)^(1 / (p - 1))
  return (y)
}

### Dual-power distortion risk measure - p in (1, infty)
# distortion function
gDPower <- function(x, p = 2) {
  y <- 1 - (1 - x)^p
  y[x >= 1] <- 1
  y[x <= 0] <- 0
  return (y)
}

# derivative of distortion function
DgDPower <- function(x, p = 2) {
  y <- p * (1 - x)^(p - 1)
  y[x > 1] <- 0
  y[x < 0] <- 0
  return (y)
}

# inverse of derivative of distortion function
DgDPowerInv <- function(x, p = 2) {
  y <- rep(0, length(x))
  ind_dom <- which((x >= 0) & (x <= p))
  y[ind_dom] <- 1 - (x[ind_dom] / p)^(1 / (p - 1))
  return (y)
}

### Wang distortion risk measure, p in (0.5, 1)
# distortion function
gWang <- function(x, p = 0.75) {
  y <- rep(0, length(x))
  ind_dom <-  which((x >= 0) & (x <= 1))
  y[ind_dom] <- pnorm(qnorm(x[ind_dom]) + qnorm(p))
  y[x > 1] <- 1
  return (y)
}

# derivative of distortion function
DgWang <- function(x, p = 0.75) {
  y <- rep(0, length(x))
  ind_dom <-  which((x >= 0) & (x <= 1))
  y[ind_dom] <- exp(-qnorm(x) * qnorm(p) - 0.5 * qnorm(p)^2)
  return (y)
}

# inverse of derivative of distortion function
DgWangInv <- function(x, p = 0.75) {
  y <- rep(0, length(x))
  ind_dom <- which(x >= 0)
  y[ind_dom] <- pnorm(-(0.5 * qnorm(p)^2 + log(x[ind_dom])) / (qnorm(p)))
  return (y)
}



###############################################
### General solution: first and k-th moment ###
###############################################

### computes distribution function given the eta vector and (g')^(-1)
# in the paper this is F^star for given eta; hence in the code F_eta is a distribution function
F_eta <- function(x, eta, DrhoInv, k = 2) {
  # x         : input values (in (0, 1))
  # eta       : vector of two eta-values
  # DrhoInv   : inverse of the derivative of the distortion function at hand
  # k         : indicating that the first and k-th moment are to be used
  
  x <- as.vector(x)                                                             # cast x into vector form
  nx <- length(x)                                                               # number of x values
  y <- rep(0, nx)                                                               # initialise functionvalues
  
  xind <- order(x)                                                              # order x
  xsort <- x[xind]                                                              # sort x
  ysort <- y                                                                    # initialise sorted y
  
  ysort[xsort < 0] <- 0                                                         # x < 0
  ysort[xsort >= 1] <- 1                                                        # x >= 1
  
  ind_dom <- which((xsort >= 0) & (xsort < 1))                                 ## x in [0, 1)
  if (length(ind_dom) > 0) {
    G_eta <- function(x) 1 - DrhoInv(eta[1] + eta[2] * x^(k - 1))               # define function G_eta, corresponding with F_eta on most places
    
    if (eta[2] <= 0) {                                                          # F_eta is constant
      ysort[ind_dom] <- 0
    } else {                                                                    # F_eta not constant
      y01 <- G_eta(xsort[ind_dom])                                              # function values of x in [0, 1)
      y01[y01 < 0] <- 0                                                         # floor at 0
      y01[y01 > 1] <- 1                                                         # ceil at 1
      indmin <- which.min(y01)                                                  # index for minimum y value
      if (y01[1] == 1) y01[1:(indmin - 1)] <- 0                                 # set first 1s to 0 if necessary
      ysort[ind_dom] <- y01                                                     # fill back in in full y
    }
  }
  
  # unsort ysort to y
  y[xind] <- ysort
  
  return (y)
}


### computes mean and k-th moment given eta vector - numerical approximation
eta2moments <- function(eta, DrhoInv, k = 2, ...) {
  # eta       : vector of two eta-values
  # DrhoInv   : inverse of the derivative of the distortion function at hand
  # k         : indicating that the first and k-th moment are to be used
  
  if (hasArg(tol)) tol <- list(...)$tol else tol = 1e-12
  if (hasArg(eps)) eps <- list(...)$eps else eps = 1e-05
  
  x <- seq(0, 1, eps)
  y <- F_eta(x, eta, DrhoInv, k)
  n <- length(x)
  
  xflatInd <- which((y[2:n] - y[1:(n - 1)]) == 0)
  yflatUnique <- unique(y[xflatInd])
  
  moms <- rep(0, 2)
  
  if (length(yflatUnique) > 0) {
    # integrate the flat parts
    for (ii in 1:length(yflatUnique)) {
      xnowind <- which(y == yflatUnique[ii])
      moms[1] = moms[1] + yflatUnique[ii] * (x[xnowind[length(xnowind)]] - x[xnowind[1]])
      moms[2] = moms[2] + yflatUnique[ii] * (x[xnowind[length(xnowind)]]^k - x[xnowind[1]]^k) / k
    }
    
    # integrate non-flat part
    if (length(x[-xflatInd]) > 2) {
      lb <- min(x[-xflatInd])
      ub <- max(x[-c(xflatInd, n)]) + eps
      moms[1] = moms[1] + integrate(f = function(x) F_eta(x, eta, DrhoInv, k), 
                                    subdivisions = 1000L, lower = lb, upper = ub,
                                    rel.tol = tol, abs.tol = tol)$value
      moms[2] = moms[2] + integrate(f = function(x) x^(k - 1) * F_eta(x, eta, DrhoInv, k), 
                                    subdivisions = 1000L, lower = lb, upper = ub,
                                    rel.tol = tol, abs.tol = tol)$value
    }
  } else {
    # integrate non-flat part (full function)
    lb <- 0
    ub <- 1
    moms[1] = moms[1] + integrate(f = function(x) F_eta(x, eta, DrhoInv, k), 
                                  subdivisions = 1000L, lower = lb, upper = ub,
                                  rel.tol = tol, abs.tol = tol)$value
    moms[2] = moms[2] + integrate(f = function(x) x^(k - 1) * F_eta(x, eta, DrhoInv, k), 
                                  subdivisions = 1000L, lower = lb, upper = ub,
                                  rel.tol = tol, abs.tol = tol)$value
  }
  
  # convert to moments
  moms <- 1 - c(1, k) * moms
  
  return (moms)
}


### compute eta values from the mean and k-th moment
moments2eta <- function(moms, DrhoInv, k = 2, etastart = c(0, 1), ...) {
  # moms      : vector of 1st and kth moment
  # DrhoInv   : inverse of the derivative of the distortion function at hand
  # k         : indicating that the first and k-th moment are to be used
  # etastart  : starting value for the optimization (using nleqslv)
  
  if (hasArg(tol)) tol <- list(...)$tol else tol = 1e-12
  if (hasArg(eps)) eps <- list(...)$eps else eps = 1e-05
  
  # make function to optimize over
  objfun <- function(eta) {
    currentmoments <- eta2moments(eta, DrhoInv, k, eps = eps, tol = tol)
    return (currentmoments - moms)
  }
  
  # solve for eta
  sol <- nleqslv(x = etastart, fn = objfun, method = "Newton",
                 control = list("xtol" = tol, "ftol" = tol))
  eta <- sol$x
  
  return (eta)
}


########################################################
### Evaluate distribution in distortion risk measure ###
########################################################

evalDistortion <- function(fDistortion, fDistribution, ...) {
  # fDistortion   : distortion function
  # fDistribution : cumulative distribution function on [0, 1]
  
  if (hasArg(tol)) tol <- list(...)$tol else tol = 1e-12
  if (hasArg(eps)) eps <- list(...)$eps else eps = 1e-05
  
  x <- seq(0, 1, eps)
  y <- fDistortion(1 - fDistribution(x))
  n <- length(x)
  
  xflatInd <- which((y[2:n] - y[1:(n - 1)]) == 0)
  yflatUnique <- unique(y[xflatInd])
  
  val <- 0
  
  if (length(yflatUnique) > 0) {
    # integrate the flat parts
    for (ii in 1:length(yflatUnique)) {
      xnowind <- which(y == yflatUnique[ii])
      val = val + yflatUnique[ii] * (x[xnowind[length(xnowind)]] - x[xnowind[1]])
    }
    
    # integrate non-flat part
    if (length(x[-xflatInd]) > 2) {
      lb <- min(x[-xflatInd])
      ub <- max(x[-c(xflatInd, n)]) + eps
      val = val + integrate(f = function(x) fDistortion(1 - fDistribution(x)), 
                            subdivisions = 1000L, lower = lb, upper = ub,
                            rel.tol = tol, abs.tol = tol)$value
    }
  } else {
    # integrate non-flat part (full function)
    lb <- 0
    ub <- 1
    val = val + integrate(f = function(x) fDistortion(1 - fDistribution(x)), 
                          subdivisions = 1000L, lower = lb, upper = ub,
                          rel.tol = tol, abs.tol = tol)$value
  }
  
  return (val)
}


#######################
### Compute moments ###
#######################

### compute m-th moment from a distribution function
distribution2moment <- function(fDistribution, m = 1, ...) {
  # fDistribution : cumulative distribution function on [0, 1]
  # m             : which moment to compute
  
  if (hasArg(tol)) tol <- list(...)$tol else tol = 1e-12
  if (hasArg(eps)) eps <- list(...)$eps else eps = 1e-05
  
  x <- seq(0, 1, eps)
  y <- fDistribution(x)
  n <- length(x)
  
  xflatInd <- which((y[2:n] - y[1:(n - 1)]) == 0)
  yflatUnique <- unique(y[xflatInd])
  
  mom <- 0
  
  if (length(yflatUnique) > 0) {
    # integrate the flat parts
    for (ii in 1:length(yflatUnique)) {
      xnowind <- which(y == yflatUnique[ii])
      mom = mom + yflatUnique[ii] * (x[xnowind[length(xnowind)]]^m - x[xnowind[1]]^m) / m
    }
    
    # integrate non-flat part
    if (length(x[-xflatInd]) > 2) {
      lb <- min(x[-xflatInd])
      ub <- max(x[-c(xflatInd, n)]) + eps
      mom = mom + integrate(f = function(x) x^(m - 1) * fDistribution(x), 
                            subdivisions = 1000L, lower = lb, upper = ub,
                            rel.tol = tol, abs.tol = tol)$value
    }
  } else {
    # integrate non-flat part (full function)
    lb <- 0
    ub <- 1
    mom = mom + integrate(f = function(x) x^(m - 1) * fDistribution(x), 
                          subdivisions = 1000L, lower = lb, upper = ub,
                          rel.tol = tol, abs.tol = tol)$value
  }
  
  # convert to moments
  mom <- 1 - m * mom
  
  return (mom)
}


###############################################################
### General solution: unbounded support - mean and variance ###
###############################################################

get_worst_case_unbounded <- function(mu, s, Drho, DrhoInv, ...) {
  
  if (hasArg(tol)) tol <- list(...)$tol else tol = 1e-12
  if (hasArg(eps)) eps <- list(...)$eps else eps = 0

  # integrate squared derivative of distortion function
  fn <- function(x) Drho(1 - x)^2
  int_g <- integrate(f = fn, lower = 0, upper = 1 - eps, subdivisions = 1000L,
                     rel.tol = tol, abs.tol = tol)$value
  
  # compute constants
  eta2 <- sqrt(int_g - 1) / s
  eta1 <- 1 - eta2 * mu
  
  # worst case distortion
  w_g <- sqrt(int_g - 1)
  Hg <- mu + w_g * s
  
  # worst case distribution function
  Fstar <- function(x) {
    
    y <- 1 - DrhoInv(eta1 + eta2 * x)
    ind_min_diff <- which.min(diff(y))
    if (y[ind_min_diff] > 0.99) y[1:ind_min_diff] <- 0
    
    y[y > 1] <- 1
    y[y < 0] <- 0
    
    return (y)
  }
  
  return (list("Hg" = Hg, "Fstar" = Fstar))
}




