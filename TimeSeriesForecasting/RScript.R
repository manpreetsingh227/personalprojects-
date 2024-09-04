# Required packages
suppressPackageStartupMessages({
  library(VineCopula) # Copula analysis
  library(ADGofTest) # Anderson-Darling Goodness-of-Fit test
  library(KScorrect) # (Lilliefors-Corrected) Kolmogorov-Smirnov Goodness-of-Fit test
  library(fGarch) # Time series analysis
  library(tseries) # Time series analysis
  library(zoo) # Data preparation
})

# Get indices prices from yahoo
FTSE100 <- get.hist.quote(instrument = "^FTSE", start = "2000-01-01", end = "2018-12-31", quote = "Adjusted", provider = "yahoo", compression = "w")
SP500 <- get.hist.quote(instrument = "^GSPC", start = "2000-01-01", end = "2018-12-31", quote = "Adjusted", provider = "yahoo", compression = "w")
SSE <- get.hist.quote(instrument = "000001.SS", start = "2000-01-01", end = "2018-12-31", quote = "Adjusted", provider = "yahoo", compression = "w")
DAX <- get.hist.quote(instrument = "^GDAXI", start = "2000-01-01", end = "2018-12-31", quote = "Adjusted", provider = "yahoo", compression = "w")
Nikkei225 <- get.hist.quote(instrument = "^N225", start = "2000-01-01", end = "2018-12-31", quote = "Adjusted", provider = "yahoo", compression = "w")
CAC40 <- get.hist.quote(instrument = "^FCHI", start = "2000-01-01", end = "2018-12-31", quote = "Adjusted", provider = "yahoo", compression = "w")

# Return names of indices with no missing values
names(Filter(function(x)all(!is.na(x)), list(FTSE100 = FTSE100, SP500 = SP500, SSE = SSE, DAX = DAX, Nikkei225 = Nikkei225, CAC40 = CAC40)))

# Calculate FTSE100 weekly log-returns
FTSE100 <- as.data.frame(coredata(diff(log(FTSE100))))
colnames(FTSE100) <- "FTSE100"

# Calculate CAC40 weekly log-returns
CAC40 <- as.data.frame(coredata(diff(log(CAC40))))
colnames(CAC40) <- "CAC40"

# Combine the log-returns data into a single dataset
Indices_log_returns <- cbind(FTSE100, CAC40)
save(Indices_log_returns, file = "Indices_log_returns.RData")

# Task 2 (a)
FTSE100 <- Indices_log_returns$FTSE100
CAC40 <- Indices_log_returns$CAC40

# Calculate mean, standard deviation, and correlation
FTSE100_mu <- mean(FTSE100)
FTSE100_sigma <- sd(FTSE100)
CAC40_mu <- mean(CAC40)
CAC40_sigma <- sd(CAC40)
correlation <- cor(FTSE100, CAC40)

# Calculate the portfolio mean and variance
port_mu <- 1/2 * (FTSE100_mu + CAC40_mu)
port_variance <- (1/2)^2 * (FTSE100_sigma^2 + CAC40_sigma^2) + 2 * 1/2 * 1/2 * correlation * FTSE100_sigma * CAC40_sigma

# Estimated 99% and 95% VaR
port_VaR <- qnorm(c(0.01, 0.05), mean = port_mu, sd = sqrt(port_variance))
message("Estimated 99% 1-week VaR: ", port_VaR[1], "\n", "Estimated 95% 1-week VaR: ", port_VaR[2])

# Scatter plots after Task 2 (a)
suppressPackageStartupMessages({library(MASS)})
set.seed(888)
mu <- c(FTSE100_mu, CAC40_mu)
sigma <- cov(Indices_log_returns)
data <- mvrnorm(n = 991, mu = mu, Sigma = sigma)
par(mfrow = c(2, 1))
plot(data, main = "Scatter plot of simulated Bivariate Normal Distribution using portfolio statistics", xlim = c(-0.3, 0.2), ylim = c(-0.3, 0.2), xlab = "We simulate a bivariate normal distribution under parametric approach's assumption.", ylab = "")
plot(FTSE100, CAC40, main = "Scatter plot of real FTSE100 and CAC40", xlab = "We observe more dispersing tails and more extreme values in our data.", ylab = "", xlim = c(-0.3, 0.2), ylim = c(-0.3, 0.2))
par(mfrow = c(1, 1))

# Task 2 (b)
jarque.bera.test(FTSE100)
jarque.bera.test(CAC40)

# Check for AR and GARCH effects
par(mfrow = c(2, 3))
acf(FTSE100, col = "green", lwd = 2)
pacf(FTSE100, col = "green", lwd = 2)
acf(FTSE100^2, col = "red", lwd = 2)
acf(CAC40, col = "green", lwd = 2)
pacf(CAC40, col = "green", lwd = 2)
acf(CAC40^2, col = "red", lwd = 2)
par(mfrow = c(1, 1))

# For loop to determine models with the lowest AIC and BIC for FTSE100
FTSE100_ICs_results <- data.frame(Lowest_AIC_Model = c("","","","","",""), AIC = c(0, 0, 0, 0, 0, 0), Lowest_BIC_Model = c("","","","","",""), BIC = c(0, 0, 0, 0, 0, 0), row.names = c("Normal (norm)", "Skew Normal (snorm)", "Student-t (std)", "Skew Student-t (sstd)", "Generalized Error (ged)", "Skew Generalized Error (sged)"))

# Normal
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'norm')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[1, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[1, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[1, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[1, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Normal
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'snorm')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[2, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[2, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[2, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[2, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Student-t
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'std')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[3, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[3, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[3, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[3, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Student-t
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'sstd')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[4, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[4, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[4, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[4, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Generalized Error
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'ged')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[5, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[5, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[5, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[5, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Generalized Error
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'sged')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[6, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[6, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[6, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[6, 4] <- ICs[which.min(ICs$BIC),][[3]]

FTSE100_ICs_results

# For loop to determine models with the lowest AIC and BIC for CAC40
CAC40_ICs_results <- data.frame(Lowest_AIC_Model = c("","","","","",""), AIC = c(0, 0, 0, 0, 0, 0), Lowest_BIC_Model = c("","","","","",""), BIC = c(0, 0, 0, 0, 0, 0), row.names = c("Normal (norm)", "Skew Normal (snorm)", "Student-t (std)", "Skew Student-t (sstd)", "Generalized Error (ged)", "Skew Generalized Error (sged)"))

# Normal
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = CAC40, trace = F, cond.dist = 'norm')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
CAC40_ICs_results[1, 1] <- ICs[which.min(ICs$AIC),][[1]]
CAC40_ICs_results[1, 2] <- ICs[which.min(ICs$AIC),][[2]]
CAC40_ICs_results[1, 3] <- ICs[which.min(ICs$BIC),][[1]]
CAC40_ICs_results[1, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Normal
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = CAC40, trace = F, cond.dist = 'snorm')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
CAC40_ICs_results[2, 1] <- ICs[which.min(ICs$AIC),][[1]]
CAC40_ICs_results[2, 2] <- ICs[which.min(ICs$AIC),][[2]]
CAC40_ICs_results[2, 3] <- ICs[which.min(ICs$BIC),][[1]]
CAC40_ICs_results[2, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Student-t
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = CAC40, trace = F, cond.dist = 'std')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
CAC40_ICs_results[3, 1] <- ICs[which.min(ICs$AIC),][[1]]
CAC40_ICs_results[3, 2] <- ICs[which.min(ICs$AIC),][[2]]
CAC40_ICs_results[3, 3] <- ICs[which.min(ICs$BIC),][[1]]
CAC40_ICs_results[3, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Student-t
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = CAC40, trace = F, cond.dist = 'sstd')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
CAC40_ICs_results[4, 1] <- ICs[which.min(ICs$AIC),][[1]]
CAC40_ICs_results[4, 2] <- ICs[which.min(ICs$AIC),][[2]]
CAC40_ICs_results[4, 3] <- ICs[which.min(ICs$BIC),][[1]]
CAC40_ICs_results[4, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Generalized Error
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = CAC40, trace = F, cond.dist = 'ged')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
CAC40_ICs_results[5, 1] <- ICs[which.min(ICs$AIC),][[1]]
CAC40_ICs_results[5, 2] <- ICs[which.min(ICs$AIC),][[2]]
CAC40_ICs_results[5, 3] <- ICs[which.min(ICs$BIC),][[1]]
CAC40_ICs_results[5, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Generalized Error
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = CAC40, trace = F, cond.dist = 'sged')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
CAC40_ICs_results[6, 1] <- ICs[which.min(ICs$AIC),][[1]]
CAC40_ICs_results[6, 2] <- ICs[which.min(ICs$AIC),][[2]]
CAC40_ICs_results[6, 3] <- ICs[which.min(ICs$BIC),][[1]]
CAC40_ICs_results[6, 4] <- ICs[which.min(ICs$BIC),][[3]]

CAC40_ICs_results

# FTSE100 model
model1 <- garchFit(formula = ~ arma(3, 0) + garch(1, 1), data = FTSE100, trace = F, cond.dist = "sstd")

# CAC40 model
model2 <- garchFit(formula = ~ arma(1, 0) + garch(1, 1), data = CAC40, trace = F, cond.dist = "sstd")

# Residuals
res1 <- residuals(model1, standardize = TRUE)
res2 <- residuals(model2, standardize = TRUE)

# Plots to check for elimination of autocorrelation and GARCH effects
par(mfrow = c(2, 2))
acf(res1, col = "green", lwd = 2)
acf(res1^2, col = "red", lwd = 2)
acf(res2, col = "green", lwd = 2)
acf(res2^2, col = "red", lwd = 2)
par(mfrow = c(1, 1))

# Box-Ljung Test
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 3)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 3)
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)

# Do the for loop again
FTSE100_ICs_results <- data.frame(Lowest_AIC_Model = c("","","","","",""), AIC = c(0, 0, 0, 0, 0, 0), Lowest_BIC_Model = c("","","","","",""), BIC = c(0, 0, 0, 0, 0, 0), row.names = c("Normal (norm)", "Skew Normal (snorm)", "Student-t (std)", "Skew Student-t (sstd)", "Generalized Error (ged)", "Skew Generalized Error (sged)"))

# Normal
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'norm')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[1, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[1, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[1, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[1, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Normal
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'snorm')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[2, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[2, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[2, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[2, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Student-t
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'std')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[3, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[3, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[3, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[3, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Student-t
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(1, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'sstd')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[4, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[4, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[4, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[4, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Generalized Error
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'ged')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[5, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[5, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[5, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[5, 4] <- ICs[which.min(ICs$BIC),][[3]]

# Skew Generalized Error
ICs <- data.frame(Model = "", AIC = 0, BIC = 0)
q <- 0
for (k in c(3, 4, 7)){
  for (i in 1:3){
    for (j in 1:3){
      q <- q + 1
      fit <- garchFit(substitute(~ arma(r, 0) + garch(p, q), list(p = i, q = j, r = k)), data = FTSE100, trace = F, cond.dist = 'sged')
      ICs <- rbind(ICs, data.frame(Model = paste("AR(", k, ")", "+ garch(", i, ",", j, ")"), AIC = fit@fit$ics[[1]], BIC = fit@fit$ics[[2]]))
    }
  }
}
FTSE100_ICs_results[6, 1] <- ICs[which.min(ICs$AIC),][[1]]
FTSE100_ICs_results[6, 2] <- ICs[which.min(ICs$AIC),][[2]]
FTSE100_ICs_results[6, 3] <- ICs[which.min(ICs$BIC),][[1]]
FTSE100_ICs_results[6, 4] <- ICs[which.min(ICs$BIC),][[3]]

FTSE100_ICs_results

# Updated model for FTSE100 and test for autocorrelation
model1 <- garchFit(formula = ~ arma(4, 0) + garch(1, 1), data = FTSE100, trace = F, cond.dist = "sstd")
res1 <- residuals(model1, standardize = TRUE)
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 4)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 4)

# Transform res1 and res2
u1 <- psstd(res1, mean = 0, sd = 1, nu = model1@fit$par["shape"], xi = model1@fit$par["skew"])[5:length(FTSE100)]
u2 <- psstd(res2, mean = 0, sd = 1, nu = model2@fit$par["shape"], xi = model2@fit$par["skew"])[5:length(CAC40)]

# Histograms
par(mfrow = c(2, 1))
hist(u1)
hist(u2)
par(mfrow = c(1, 1))

# Kolmogorov-Smirnov test
KStest1 <- LcKS(u1, cdf = "punif")
KStest1$p.value
KStest2 <- LcKS(u2, cdf = "punif")
KStest2$p.value

# Anderson-Darling test
ADtest1 <- ad.test(u1, null = "punif")
ADtest1$p.value
ADtest2 <- ad.test(u2, null = "punif")
ADtest2$p.value

N <- 10000
set.seed(888)
u_sim <- BiCopSim(N, obj = model)
# Apply Inverse Probability Integral Transform (IPIT) to the simulated results for both indices.
res1_sim <- qsstd(u_sim[,1], mean = 0, sd = 1, nu = model1@fit$par["shape"], xi = model1@fit$par["skew"])
res2_sim <- qsstd(u_sim[,2], mean = 0, sd = 1, nu = model2@fit$par["shape"], xi = model2@fit$par["skew"])

# Compare the original and simulated residuals
par(mfrow = c(2, 1))
plot(res1, res2, xlim = c(-9, 5), ylim = c(-11, 6), main = "Scatter plot of residuals of Model 1 and Model 2", xlab = "Residuals of Model 1", ylab = "Residuals of Model 2")
plot(res1_sim, res2_sim, xlim = c(-9, 5), ylim = c(-11, 6), main = "Scatter plot of simulated residuals", xlab = "Simulated residuals 1", ylab = "Simulated residuals 2")
par(mfrow = c(1, 1))

## FTSE100
# Extract parameters from model1 for later calculations
mu_1 <- model1@fit$par["mu"]
ar1_1 <- model1@fit$par["ar1"]
ar2_1 <- model1@fit$par["ar2"]
ar3_1 <- model1@fit$par["ar3"]
ar4_1 <- model1@fit$par["ar4"]
omega_1 <- model1@fit$par["omega"]
alpha_1 <- model1@fit$par["alpha1"]
beta_1 <- model1@fit$par["beta1"]

# Create a vector to store the simulated 1-week ahead log-returns
y1_sim <- numeric(N)

# Reintroduce the GARCH effects
sigma2s_1 <- omega_1 + alpha_1 * model1@residuals[991]^2 + beta_1 * model1@sigma.t[991]^2

# Reintroduce the AR effects and store the simulated results
for (i in 1:N) {
  y1_sim[i] <- mu_1 + ar1_1 * FTSE100[991] + ar2_1 * FTSE100[990] + ar3_1 * FTSE100[989] + ar4_1 * FTSE100[988] + sqrt(sigma2s_1) * res1_sim[i]
}

## CAC40
# Extract parameters from model2 for later calculations
mu_2 <- model2@fit$par["mu"]
ar1_2 <- model2@fit$par["ar1"]
omega_2 <- model2@fit$par["omega"]
alpha_2 <- model2@fit$par["alpha1"]
beta_2 <- model2@fit$par["beta1"]

# Create a vector to store the simulated 1-week ahead log-returns
y2_sim <- numeric(N)

# Reintroduce the GARCH effects
sigma2s_2 <- omega_2 + alpha_2 * model2@residuals[991]^2 + beta_2 * model2@sigma.t[991]^2

# Reintroduce the AR effects and store the simulated results
for (i in 1:N) {
  y2_sim[i] <- mu_2 + ar1_2 * CAC40[991] + sqrt(sigma2s_2) * res2_sim[i]
}

# Calculate portfolio VaR
port_sim <- matrix(0, nrow = N, ncol = 1)
VaR_sim <- matrix(0, nrow = 1, ncol = 2)

port_sim <- log(1 + ((exp(y1_sim) - 1) + (exp(y2_sim) - 1)) * (1/2))
VaR_sim <- quantile(port_sim, c(0.01, 0.05))
message("Estimated 99% 1-week VaR: ", VaR_sim[1], "\n", "Estimated 95% 1-week VaR: ", VaR_sim[2])
