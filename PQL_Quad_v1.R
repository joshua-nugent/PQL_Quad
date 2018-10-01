library(lme4)
library(MASS)
library(tidyverse)

check_estimates <- function(n, p, b_0, b_1, sigbsq){
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
  
  glmmresPQL <- glmmPQL(y ~ x1, random = ~ 1 | id , family = binomial(link = "logit"))
  glmmresLP <- glmer(y ~ x1 + (1|id), nAGQ = 1, family = binomial(link = "logit"))
  glmmres4 <- glmer(y ~ x1 + (1|id), nAGQ = 4, family = binomial(link = "logit"))
  glmmres10 <- glmer(y ~ x1 + (1|id), nAGQ = 10, family = binomial(link = "logit"))
  glmmres25 <- glmer(y ~ x1 + (1|id), nAGQ = 25, family = binomial(link = "logit"))
  
  betaPQL <- fixef(glmmresPQL)
  betaLP <- fixef(glmmresLP)
  beta4 <- fixef(glmmres4)
  beta10 <- fixef(glmmres10)
  beta25 <- fixef(glmmres25)
 
  value_labels <- c("n", "p", "b_0", "b_1", "sigbsq", 
              "PQL_b_0", "Laplace_b_0", "nAGQ_4_b_0", "nAGQ_10_b_0", "nAGQ_25_b_0",
              "PQL_b_1", "Laplace_b_1", "nAGQ_4_b_1", "nAGQ_10_b_1", "nAGQ_25_b_1")
  estimates <-c(n, p, b_0, b_1, sigbsq,
                betaPQL[1], betaLP[1], beta4[1], beta10[1], beta25[1],
                betaPQL[2], betaLP[2], beta4[2], beta10[2], beta25[2])

  names(estimates) <- value_labels
  return(estimates)
  
####  Might use this later for later SE / coverage analysis 
#  VcovPQL <- vcov(glmmresPQL, useScale = FALSE) - get covariance matrix...
# ...and use diagonal entries to get SE
#  sePQL <- sqrt(diag(VcovPQL))
}


###############################################################
# Here we actually start running the simulation
###############################################################

# Parameters to input to the generator...
# We will explore what gives interesting results.
n <- 40
p <- 100
b_0 <- -2
b_1 <- 0.5
sigbsq <- 1

# KEY STEP - Output of the models
data <- as_tibble(t(replicate(1000, check_estimates(n, p, b_0, b_1, sigbsq))))

# Add columns to our generator output to measure bias of b_0 term...
# ...later, this could be put into the function as well.
data <- data %>% mutate(
  b_0_PQL_bias = b_0 - PQL_b_0,
  b_0_Laplace_bias = b_0 - Laplace_b_0,
  b_0_4_bias = b_0 - nAGQ_4_b_0,
  b_0_10_bias = b_0 - nAGQ_10_b_0,
  b_0_25_bias = b_0 - nAGQ_25_b_0
  )
# Add columns to our generator output to measure bias of b_1 term...
data <- data %>% mutate(
  b_1_PQL_bias = b_1 - PQL_b_1,
  b_1_Laplace_bias = b_1 - Laplace_b_1,
  b_1_4_bias = b_1 - nAGQ_4_b_1,
  b_1_10_bias = b_1 - nAGQ_10_b_1,
  b_1_25_bias = b_1 - nAGQ_25_b_1
)

# Organizational structure as model list might grow
# to hopefully avoid repeated code
diff_col_names_b_0 <- c(
  "b_0_PQL_bias",
  "b_0_Laplace_bias",
  "b_0_4_bias",
  "b_0_10_bias",
  "b_0_25_bias"
  )
diff_col_names_b_1 <- c(
  "b_1_PQL_bias",
  "b_1_Laplace_bias",
  "b_1_4_bias",
  "b_1_10_bias",
  "b_1_25_bias"
)


# Look at bias / variance
for (i in diff_col_names_b_0){
  cat(sprintf("mean for intercept error of %s", i))
  cat(sprintf("is: \n"))
  m <- mean(data[[i]])
  print(m)
}
for (i in diff_col_names_b_0){
  cat(sprintf("std dev for intercept error of %s", i))
  cat(sprintf("is: \n"))
  sd <- sd(data[[i]])
  print(sd)
}
for (i in diff_col_names_b_0){
  histogram <- hist(data[[i]], main = paste("Histogram of intercept error of", i))
}


# Same as above, but for beta1 term
for (i in diff_col_names_b_1){
  cat(sprintf("mean for b_1 estimate error for %s", i))
  cat(sprintf("is: \n"))
  m <- mean(data[[i]])
  print(m)
}
for (i in diff_col_names_b_1){
  cat(sprintf("std dev for b_1 estimate error for %s", i))
  cat(sprintf("is: \n"))
  sd <- sd(data[[i]])
  print(sd)
}
for (i in diff_col_names_b_1){
  histogram <- hist(data[[i]], main = paste("Histogram of b_1 estimate error for", i))
}

