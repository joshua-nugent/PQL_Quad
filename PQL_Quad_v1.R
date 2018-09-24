library(lme4)
library(MASS)

# old: beta <- c(-2, 1.5, 0.5, -1)

check_estimates <- function(n, p, b_0, b_1, sigbsq){
  ## Initial tests: n <- 1500 p <- 3 b_0 <- -2 b_1 <- 1.5 sigbsq <- 4
  
  betas <- c(b_0, b_1)

  # vector of p copies of each ID for 'long' data
  id <- rep(1:n, each = p)
  
  # vector of 1111...0000 for treatment arms
  x1 <- as.numeric(id < (n+1)/2)
  
  # Generate random normal values for error term, replicated p times for each ID
  randint <- rep(rnorm(n, 0, sqrt(sigbsq)), each = p)
  
  # Data generation with specified model plus noise
  linpred <- b_0 + b_1*x1 + randint
  
  # Invert t
  expit <- exp(linpred) / (1 + exp(linpred))
  
  y <- runif(p*n) < expit
  
  # ? glmmML not needed anymore?
  glmmresPQL <- glmmPQL(y ~ x1, random = ~ 1 | id , family = binomial(link = "logit"))
  glmmresLP <- glmer(y ~ x1 + (1|id), nAGQ = 1, family = binomial(link = "logit"))
  glmmres10 <- glmer(y ~ x1 + (1|id), nAGQ = 10, family = binomial(link = "logit"))
  glmmres50 <- glmer(y ~ x1 + (1|id), nAGQ = 50, family = binomial(link = "logit"))
  glmmres100 <- glmer(y ~ x1 + (1|id), nAGQ = 100, family = binomial(link = "logit"))
  
  (betaPQL <- fixef(glmmresPQL))
  (betaLP <- fixef(glmmresLP))
  (beta10 <- fixef(glmmres10))
  (beta50 <- fixef(glmmres50))
  (beta100 <- fixef(glmmres100))
 
  method <- c("True Value", "PQL", "Laplace", "nAGQ_10", "nAGQ_50", "nAGQ_100")
  intercept_estimate <-c(b_0, betaPQL[1], betaLP[1], beta10[1], beta50[1], beta100[1])
  x1_estimate <-c(b_1, betaPQL[2], betaLP[2], beta10[2], beta50[2], beta100[2])

  results <- data.frame(method, intercept_estimate, x1_estimate)
  names(results) <- c("Method", "Intercept Estimate", "Beta1 Estimate")
  
  return(results)
  
####  For later SE / coverage analysis 
#  VcovPQL <- vcov(glmmresPQL, useScale = FALSE)
#  VcovLP <- vcov(glmmresLP, useScale = FALSE)
#  Vcov10 <- vcov(glmmres10, useScale = FALSE)
#  Vcov50 <- vcov(glmmres50, useScale = FALSE)
#  Vcov100 <- vcov(glmmres100, useScale = FALSE)
  
#  (sePQL <- sqrt(diag(VcovPQL)))
#  (seLP <- sqrt(diag(VcovLP)))
#  (se10 <- sqrt(diag(Vcov10)))
#  (se50 <- sqrt(diag(Vcov50)))
#  (se100 <- sqrt(diag(Vcov100)))
}

check_estimates(1500, 3, -2, 1.5, 4)
check_estimates(150, 3, -2, 1.5, 4)
check_estimates(1500, 3, -2, 1.5, 1)
check_estimates(1500, 3, -2, 1.5, 10)
