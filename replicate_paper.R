rm(list = ls())


### Make sure to set your working directory to the one containing this file and "functions_paper.R"
source('functions_paper.R')


##################
### Remark 2.4 ###
##################

### Define distortion function
p <- 2                                                      # set exponent
dfun <- function(x) gDPower(x, p = p)                       # define distortion function
DdfunInv <- function(x) DgDPowerInv(x, p = p)               # define inverse of the derivative of the distortion function

### plot distortion function
x <- seq(0, 1, 0.001)
plot(x, dfun(x), type = 'l', lwd = 3,
     xlab = "x", ylab = "g(x)", cex.lab = 1.7, cex.axis = 1.7, mgp = c(2.7, 1, 0))

### inequalities counterexample
moms <- c(0.5, 0.33, 0.24)                                  # first three moments

eta12 <- moments2eta(moms[c(1, 2)], DdfunInv, 2)            # find eta_12 when fixing 1st and 2nd moment
eta13 <- moments2eta(moms[c(1, 3)], DdfunInv, 3)            # find eta_13 when fixing 1st and 3th moment

moms12 <- moms13 <- rep(0, 3)                               # initialize moments of the solution cdf's
moms12[1] <- distribution2moment(function(x) F_eta(x, eta12, DdfunInv, 2), 1) # 1st moment of F_eta12
moms12[2] <- distribution2moment(function(x) F_eta(x, eta12, DdfunInv, 2), 2) # 2nd moment of F_eta12
moms12[3] <- distribution2moment(function(x) F_eta(x, eta12, DdfunInv, 2), 3) # 3rd moment of F_eta12
moms13[1] <- distribution2moment(function(x) F_eta(x, eta13, DdfunInv, 3), 1) # 1st moment of F_eta13
moms13[2] <- distribution2moment(function(x) F_eta(x, eta13, DdfunInv, 3), 2) # 2nd moment of F_eta13
moms13[3] <- distribution2moment(function(x) F_eta(x, eta13, DdfunInv, 3), 3) # 3rd moment of F_eta13

round(moms12, 4)
round(moms13, 4)

plot(x, F_eta(x, eta12, DdfunInv, 2), type = 'l', lwd = 3, lty = 2,
     col = "black", ylim = c(0, 1), cex.lab = 1.7, cex.axis = 1.7,
     xlab = "x", ylab = "F(x)", mgp = c(2.7, 1, 0))
lines(x, F_eta(x, eta13, DdfunInv, 3), type = 'l', lwd = 3, col = "blue")
legend("topleft", c("1st and 2nd moment", "1st and 3th moment"), lwd = c(3, 3), lty = c(2, 1),
       col = c("black", "blue"), cex = 1.7)

### F* as in (14)
# F_12
eta12star <- round(eta12 / 2, 4)
round((c(0, 1) - eta12star[1]) / eta12star[2], 4)           # bounds where F_12(x) = eta12star[1] + eta12star[2] * x

# F_13
eta13star <- round(eta13 / 2, 4)
round(sqrt((c(0, 1) - eta13star[1]) / eta13star[2]), 4)     # bounds where F_13(x) = eta13star[1] + eta13star[2] * x^2

### Distortion values
evalDistortion(dfun, function(x) F_eta(x, eta12, DdfunInv, 2))
evalDistortion(dfun, function(x) F_eta(x, eta13, DdfunInv, 3))


###################
### Section 3.1 ###
###################

##### comparison of distortion functions
### define distortion functions parameters
alpha_Power <- c(1 / 2, 1 / 5, 1 / 10)
beta_DPower <- c(2, 5, 10)
p_Wang <- c(0.8, 0.9, 0.95)

### initialize variables
maxvalues <- matrix(NA, nrow = 9, ncol = 4)

eta12 <- eta13 <- eta14 <-
  matrix(NA, nrow = 9, ncol = 2)
eta123 <- matrix(NA, nrow = 9, ncol = 3)

moms12 <- moms13 <- moms14 <-
  matrix(NA, nrow = 9, ncol = 2)
moms123 <- matrix(NA, nrow = 9, ncol = 3)

moments <- c(1 / 2, 1 / 3, 1 / 4, 1 / 5)

### first and second moment fixed
k <- 2
for (ii in 1:length(alpha_Power)) {
  # Power distortion
  eta12[ii,] <- moments2eta(moments[c(1, k)], function(x) DgPowerInv(x, alpha_Power[ii]), k)
  fDist <- function(x) F_eta(x, eta12[ii,], function(x) DgPowerInv(x, alpha_Power[ii]), k)
  moms12[ii, 1] <- distribution2moment(fDist, 1)
  moms12[ii, 2] <- distribution2moment(fDist, k)
  maxvalues[ii, 1] <- evalDistortion(function(x) gPower(x, alpha_Power[ii]), fDist)
  
  maxvalues[ii, 4] <- get_worst_case_unbounded(mu = moments[1], s = sqrt(moments[2] - moments[1]^2),
                                               Drho = function(x) DgPower(x, alpha_Power[ii]),
                                               DrhoInv = function(x) DgPowerInv(x, alpha_Power[ii]), eps = 1e-6, tol = 1e-7)$Hg
  
  # Dual-power distortion
  eta12[length(alpha_Power) + ii,] <- moments2eta(moments[c(1, k)], function(x) DgDPowerInv(x, beta_DPower[ii]), k, etastart = c(-8.03e-01, 3.21), tol = 2e-07)
  fDist <- function(x) F_eta(x, eta12[length(alpha_Power) + ii,], function(x) DgDPowerInv(x, beta_DPower[ii]), k)
  moms12[length(alpha_Power) + ii, 1] <- distribution2moment(fDist, 1, tol = 1e-07)
  moms12[length(alpha_Power) + ii, 2] <- distribution2moment(fDist, k, tol = 1e-07)
  maxvalues[length(alpha_Power) + ii, 1] <- evalDistortion(function(x) gDPower(x, beta_DPower[ii]), fDist, tol = 1e-07)
  
  maxvalues[length(alpha_Power) + ii, 4] <- get_worst_case_unbounded(mu = moments[1], s = sqrt(moments[2] - moments[1]^2), 
                                                                     Drho = function(x) DgDPower(x, beta_DPower[ii]), 
                                                                     DrhoInv = function(x) DgDPowerInv(x, beta_DPower[ii]))$Hg
  
  # Wang distortion
  eta12[length(c(alpha_Power, beta_DPower)) + ii,] <- moments2eta(moments[c(1, k)], function(x) DgWangInv(x, p_Wang[ii]), k)
  fDist <- function(x) F_eta(x, eta12[length(c(alpha_Power, beta_DPower)) + ii,], function(x) DgWangInv(x, p_Wang[ii]), k)
  moms12[length(c(alpha_Power, beta_DPower)) + ii, 1] <- distribution2moment(fDist, 1)
  moms12[length(c(alpha_Power, beta_DPower)) + ii, 2] <- distribution2moment(fDist, k)
  maxvalues[length(c(alpha_Power, beta_DPower)) + ii, 1] <- evalDistortion(function(x) gWang(x, p_Wang[ii]), fDist)
  
  maxvalues[length(c(alpha_Power, beta_DPower)) + ii, 4] <- 
    get_worst_case_unbounded(mu = moments[1], s = sqrt(moments[2] - moments[1]^2), 
                             Drho = function(x) DgWang(x, p_Wang[ii]), 
                             DrhoInv = function(x) DgWangInv(x, p_Wang[ii]), tol = 1e-7)$Hg
}
print(moms12) # as check

### first and third moment fixed
k <- 3
for (ii in 1:length(alpha_Power)) {
  # Power distortion
  eta13[ii,] <- moments2eta(moments[c(1, k)], function(x) DgPowerInv(x, alpha_Power[ii]), k)
  fDist <- function(x) F_eta(x, eta13[ii,], function(x) DgPowerInv(x, alpha_Power[ii]), k)
  moms13[ii, 1] <- distribution2moment(fDist, 1)
  moms13[ii, 2] <- distribution2moment(fDist, k)
  maxvalues[ii, 2] <- evalDistortion(function(x) gPower(x, alpha_Power[ii]), fDist)
  
  # Dual-power distortion
  eta13[length(alpha_Power) + ii,] <- moments2eta(moments[c(1, k)], function(x) DgDPowerInv(x, beta_DPower[ii]), k, etastart = c(-0.20537485, 3.4658162), tol = 1e-07)
  fDist <- function(x) F_eta(x, eta13[length(alpha_Power) + ii,], function(x) DgDPowerInv(x, beta_DPower[ii]), k)
  moms13[length(alpha_Power) + ii, 1] <- distribution2moment(fDist, 1, tol = 1e-07)
  moms13[length(alpha_Power) + ii, 2] <- distribution2moment(fDist, k, tol = 1e-07)
  maxvalues[length(alpha_Power) + ii, 2] <- evalDistortion(function(x) gDPower(x, beta_DPower[ii]), fDist, tol = 1e-07)
  
  # Wang distortion
  eta13[length(c(alpha_Power, beta_DPower)) + ii,] <- moments2eta(moments[c(1, k)], function(x) DgWangInv(x, p_Wang[ii]), k)
  fDist <- function(x) F_eta(x, eta13[length(c(alpha_Power, beta_DPower)) + ii,], function(x) DgWangInv(x, p_Wang[ii]), k)
  moms13[length(c(alpha_Power, beta_DPower)) + ii, 1] <- distribution2moment(fDist, 1)
  moms13[length(c(alpha_Power, beta_DPower)) + ii, 2] <- distribution2moment(fDist, k)
  maxvalues[length(c(alpha_Power, beta_DPower)) + ii, 2] <- evalDistortion(function(x) gWang(x, p_Wang[ii]), fDist)
}
print(moms13) # as check

### first and fourth moment fixed
k <- 4
for (ii in 1:length(alpha_Power)) {
  # Power distortion
  eta14[ii,] <- moments2eta(moments[c(1, k)], function(x) DgPowerInv(x, alpha_Power[ii]), k)
  fDist <- function(x) F_eta(x, eta14[ii,], function(x) DgPowerInv(x, alpha_Power[ii]), k)
  moms14[ii, 1] <- distribution2moment(fDist, 1)
  moms14[ii, 2] <- distribution2moment(fDist, k)
  maxvalues[ii, 3] <- evalDistortion(function(x) gPower(x, alpha_Power[ii]), fDist)
  
  # Dual-power distortion
  eta14[length(alpha_Power) + ii,] <- moments2eta(moments[c(1, k)], function(x) DgDPowerInv(x, beta_DPower[ii]), k, etastart = c(-0.0349401504, 4.1730171), tol = 5e-08)
  fDist <- function(x) F_eta(x, eta14[length(alpha_Power) + ii,], function(x) DgDPowerInv(x, beta_DPower[ii]), k)
  moms14[length(alpha_Power) + ii, 1] <- distribution2moment(fDist, 1, tol = 5e-08)
  moms14[length(alpha_Power) + ii, 2] <- distribution2moment(fDist, k, tol = 5e-08)
  maxvalues[length(alpha_Power) + ii, 3] <- evalDistortion(function(x) gDPower(x, beta_DPower[ii]), fDist)
  
  # Wang distortion
  eta14[length(c(alpha_Power, beta_DPower)) + ii,] <- moments2eta(moments[c(1, k)], function(x) DgWangInv(x, p_Wang[ii]), k)
  fDist <- function(x) F_eta(x, eta14[length(c(alpha_Power, beta_DPower)) + ii,], function(x) DgWangInv(x, p_Wang[ii]), k)
  moms14[length(c(alpha_Power, beta_DPower)) + ii, 1] <- distribution2moment(fDist, 1)
  moms14[length(c(alpha_Power, beta_DPower)) + ii, 2] <- distribution2moment(fDist, k)
  maxvalues[length(c(alpha_Power, beta_DPower)) + ii, 3] <- evalDistortion(function(x) gWang(x, p_Wang[ii]), fDist)
}
print(moms14) # as check

### distortion values
round(maxvalues, 4)

### make the plots

# distortion functions
x <- seq(0, 0.9999, 0.005)
plot(x, gPower(x, 1 / 5), type = 'l', lwd = 3,
     xlab = "x", ylab = "g(x)", cex.lab = 1.7, cex.axis = 1.7, mgp = c(2.7, 1, 0))
lines(x, gDPower(x, 5), lwd = 3, col = "blue", lty = 2)
lines(x, gWang(x, 0.8), lwd = 3, col = "red", lty = 3)
legend("bottomright", c("power (alpha = 0.2)", "dual-power (beta = 5)", "Wang (p = 0.8)"), lwd = c(3, 3, 3), lty = c(1, 2, 3),
       col = c("black", "blue", "red"), cex = 1.7)

# maximizing distributions under power distortion
xf <- seq(0, 0.9999, 0.001)
plot(xf, F_eta(xf, eta12[2,], function(x) DgPowerInv(x, alpha_Power[2]), 2),
     type = 'l', lwd = 3, xlab = "x", ylab = "F(x)", cex.lab = 1.7, cex.axis = 1.7, ylim = c(0, 1), mgp = c(2.7, 1, 0))
lines(x, F_eta(x, eta13[2,], function(x) DgPowerInv(x, alpha_Power[2]), 3),
      lwd = 3, col = "blue", lty = 2)
lines(xf, F_eta(xf, eta14[2,], function(x) DgPowerInv(x, alpha_Power[2]), 4),
      lwd = 3, col = "red", lty = 3)
legend("topleft", c("1st and 2nd moment", "1st and 3rd moment", "1st and 4th moment"), 
       lwd = c(3, 3, 3), lty = c(1, 2, 3), col = c("black", "blue", "red"), cex = 1.7)

# maximizing distributions under dual-power distortion
xf <- seq(0, 0.9999, 0.001)
plot(xf, F_eta(xf, eta12[length(alpha_Power) + 2,], function(x) DgDPowerInv(x, beta_DPower[2]), 2),
     type = 'l', lwd = 3, xlab = "x", ylab = "F(x)", cex.lab = 1.7, cex.axis = 1.7, ylim = c(0, 1), mgp = c(2.7, 1, 0))
lines(x, F_eta(x, eta13[length(alpha_Power) + 2,], function(x) DgDPowerInv(x, beta_DPower[2]), 3),
      lwd = 3, col = "blue", lty = 2)
lines(xf, F_eta(xf, eta14[length(alpha_Power) + 2,], function(x) DgDPowerInv(x, beta_DPower[2]), 4),
      lwd = 3, col = "red", lty = 3)
legend("topleft", c("1st and 2nd moment", "1st and 3rd moment", "1st and 4th moment"), 
       lwd = c(3, 3, 3), lty = c(1, 2, 3), col = c("black", "blue", "red"), cex = 1.7)

# maximizing distributions under power Wang distortion
xf <- seq(0, 0.9999, 0.001)
plot(xf, F_eta(xf, eta12[length(c(alpha_Power, beta_DPower)) + 2,], function(x) DgWangInv(x, p_Wang[2]), 2),
     type = 'l', lwd = 3, xlab = "x", ylab = "F(x)", cex.lab = 1.7, cex.axis = 1.7, ylim = c(0, 1), mgp = c(2.7, 1, 0))
lines(x, F_eta(x, eta13[length(c(alpha_Power, beta_DPower)) + 2,], function(x) DgWangInv(x, p_Wang[2]), 3),
      lwd = 3, col = "blue", lty = 2)
lines(xf, F_eta(xf, eta14[length(c(alpha_Power, beta_DPower)) + 2,], function(x) DgWangInv(x, p_Wang[2]), 4),
      lwd = 3, col = "red", lty = 3)
legend("topleft", c("1st and 2nd moment", "1st and 3rd moment", "1st and 4th moment"), 
       lwd = c(3, 3, 3), lty = c(1, 2, 3), col = c("black", "blue", "red"), cex = 1.7)


###################
### Section 3.2 ###
###################

x <- seq(0, 1, 0.0001)

### Plot moment space
y.low <- x^2
y.high <- x
plot(x, y.high, type = 'n',
     xlab = expression(c[1]), ylab = expression(c[2]), cex.lab = 1.7, cex.axis = 1.7,
     cex = 1.7, las = 1, mgp = c(2.7, 1, 0))
lines(x, y.low, col = 'grey')
lines(x, y.high, col = 'grey')
polygon(c(x, rev(x)), c(y.high, rev(y.low)), col = "grey30", border = NA)

### Plot moment space - subsets
y.low <- x^2
y.high <- x
plot(x, y.high, type = 'n',
     xlab = expression(c[1]), ylab = expression(c[2]), cex.lab = 1.7, cex.axis = 1.7,
     cex = 1.7, mgp = c(2.7, 1, 0))
lines(x, y.low, col = 'grey')
lines(x, y.high, col = 'grey')
polygon(c(x, rev(x)), c(y.high, rev(y.low)), col = "grey30", border = NA)

x121 <- seq(from = 0.5, to = 1, by = 0.0001)
y.low <- (4 * x121^2 - 2 * x121 + 1) / 3
y.high <- (4 * x121 - 1) / 3
polygon(c(x121, rev(x121)), c(y.high, rev(y.low)), col = "grey80", border = NA)

x012 <- seq(from = 0, to = 0.5, by = 0.0001)
y.low <- 4 * x012^2 / 3
y.high <- 2 * x012 / 3
polygon(c(x012, rev(x012)), c(y.high, rev(y.low)), col = "grey55", border = NA)

x <- c(x012, x121[-1])
y.low <- x^2
y.high <- c(4 * x012^2 / 3, (4 * x121[-1]^2 - 2 * x121[-1] + 1) / 3)
polygon(c(x, rev(x)), c(y.high, rev(y.low)), col = "grey10", border = NA)

legend("topleft", c(expression(N[1]), expression(N[2]), expression(N[3]), expression(N[4])),
       lwd = c(10, 10, 10, 10), col = c("grey30", "grey80", "grey55", "grey10"),
       cex = 1.7)

# Plot E set - subsets
eta1 <- seq(from = 0, to = 2, by = 0.0001)
eta2 <- 2 - eta1
plot(eta1, eta2, type = 'n', xlim = c(-4, 2.5), ylim = c(-0.5, 4),
     xlab = expression(eta[1]), ylab = expression(eta[2]), cex.lab = 1.7, cex.axis = 1.7,
     cex = 1.7, mgp = c(2.7, 1, 0))
y.low <- rep(0, length(eta1))
y.high <- eta2
polygon(c(eta1, rev(eta1)), c(y.high, rev(y.low)), col = "grey30", border = NA)

eta1 <- seq(from = -4.5, to = 0, by = 0.0001)
y.low <- -eta1
y.high <- 2 - eta1
polygon(c(eta1, rev(eta1)), c(y.high, rev(y.low)), col = "grey80", border = NA)

eta1 <- seq(from = 0, to = 2, by = 0.0001)
y.low <- 2 - eta1
y.high <- rep(4.5, length(eta1))
polygon(c(eta1, rev(eta1)), c(y.high, rev(y.low)), col = "grey55", border = NA)

eta1 <- seq(from = -4.5, to = 0, by = 0.0001)
y.low <- 2 - eta1
y.high <- rep(4.5, length(eta1))
polygon(c(eta1, rev(eta1)), c(y.high, rev(y.low)), col = "grey10", border = NA)

abline(h = 0, lwd = 2, lty = 4)
abline(v = 0, lwd = 2, lty = 4)

legend(-4, 1.9, c(expression(E[1]), expression(E[2]), expression(E[3]), expression(E[4])),
       lwd = c(10, 10, 10, 10), col = c("grey30", "grey80", "grey55", "grey10"),
       cex = 1.7)

