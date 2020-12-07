### Calculating Cox-Snell R2 from McFadden's R2

CSR2_from_MFR2 <- function(R2MFapp, E, n, p){
  
  # function to recover Cox-Snell R2 from McFadden's R2
  # 
  # for use when calculating sample size for developing a prediction model (pmsampsize)
  # Riley RD, Snell KIE, Ensor S, Burke DL, Harrell Jr FE, Moons KG, Collns GS. Minimum
  # sample size for developing a multivariable prediction model: part II - binary and
  # time-to-event outcomes. Statistics in Medicine 2019; 38: 1276-1296.
  #
  # R2MF: McFadden's R2
  # E: total number of outcome events
  # n: sample size
  # p: number of predictor parameters
  
  # EQNS refer to equations reported in Riley et al (above)
  
  lnLNULL <- E * log(E/n) + (n - E) * log(1 - (E/n))  # EQN 12
  
  lnLmodel <- (1 - R2MFapp) * lnLNULL  # from EQN 16
  
  LR <- -2 * (lnLNULL - lnLmodel) # EQN 4

  R2CSapp <- 1 - exp(-LR/n) # EQN 6
  
  SVH <- 1 + (p / (n * log(1 - R2CSapp))) # EQN 7
  
  R2CSadj <- SVH * R2CSapp # EQN 8
  
  return(R2CSadj)
}

### quick simulation to see if the approximation works
require(performance)
generate_data <- function(NN, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors, cor0 = 0.1, cor1 = 0.05, beta.0 = 0){
  
  n.predictors <- n.true.predictors + n.noise.predictors
  mu0 <- rep(0, n.predictors)
  
  # Specify correlation matrix
  Sigma0 <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
  Sigma0[1:n.true.predictors, 1:n.true.predictors] <- cor0
  Sigma0[(n.true.predictors+1):n.predictors, (n.true.predictors+1):n.predictors] <- cor1
  diag(Sigma0) <- 1.0
  
  x <- mvrnorm(NN, mu0, Sigma0)
  
  beta <- c(0.5, 0.3, 0.3, 0.25, 0.25, rep(0, n.noise.predictors))
  
  y <- runif(NN) < 1 / (1 + exp(-beta.0 - x %*% beta))
  
  DATA   <- data.frame(x)
  DATA$y <- y * 1
  DATA
}

NSIM <- 100
N <- c(100, 500, 1000, 10000)
R2CSadj     <- matrix(ncol = length(N), nrow = NSIM)
R2CSadj_est <- matrix(ncol = length(N), nrow = NSIM)
for(j in 1:length(N)){
  for(i in 1:NSIM){
  # development data
    dat1 <- generate_data(NN = N[j], n.true.predictors = 5, n.noise.predictors = 1) 
    n <- nrow(dat1)
    E <- sum(dat1$y == 1)
    fit.glm <- glm(y~., data=dat1, family = 'binomial')
    tmp <- performance::r2_coxsnell(fit.glm)  # apparent Cox_snell r2 from package (as a check)
    R2CSadj[i,j] <- tmp * (1+ (length(coef(fit.glm))-1) / (n * log(1-tmp))) # need to multiply tmp by SVH to get adjusted Cox-Snell
    R2CSadj_est[i,j] <- CSR2_from_MFR2(as.numeric(r2_mcfadden(fit.glm)[1]), E = E, n = n, p = 6)
  }
}

par(mfrow = c(2, 2))
plot(R2CSadj[,1], R2CSadj_est[, 1], pch = 20, xlab = "Cox-Snell R2", ylab = "Recovered Cox-Snell R2", main="n=100 \n (50% outcome fraction)")
grid()
abline(a = 0, b = 1)

plot(R2CSadj[, 2], R2CSadj_est[, 2], pch = 20, xlab = "Cox-Snell R2", ylab = "Recovered Cox-Snell R2", main="n=500 \n (50% outcome fraction)")
grid()
abline(a = 0, b = 1)

plot(R2CSadj[, 3], R2CSadj_est[, 3], pch = 20, xlab = "Cox-Snell R2", ylab = "Recovered Cox-Snell R2", main="n=1000 \n (50% outcome fraction)")
grid()
abline(a = 0, b = 1)

plot(R2CSadj[, 4], R2CSadj_est[, 4], pch = 20, xlab = "Cox-Snell R2", ylab = "Recovered Cox-Snell R2", main="n=10000 \n (50% outcome fraction)")
grid()
abline(a = 0, b = 1)