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
  glmmres10 <- glmer(y ~ x1 + (1|id), nAGQ = 10, family = binomial(link = "logit"))
  glmmres50 <- glmer(y ~ x1 + (1|id), nAGQ = 50, family = binomial(link = "logit"))
  glmmres100 <- glmer(y ~ x1 + (1|id), nAGQ = 100, family = binomial(link = "logit"))
  
  betaPQL <- fixef(glmmresPQL)
  betaLP <- fixef(glmmresLP)
  beta10 <- fixef(glmmres10)
  beta50 <- fixef(glmmres50)
  beta100 <- fixef(glmmres100)
 
  value_labels <- c("n", "p", "b_0", "b_1", "sigbsq", 
              "PQL_b_0", "Laplace_b_0", "nAGQ_10_b_0", "nAGQ_50_b_0", "nAGQ_100_b_0",
              "PQL_b_1", "Laplace_b_1", "nAGQ_10_b_1", "nAGQ_50_b_1", "nAGQ_100_b_1")
  estimates <-c(n, p, b_0, b_1, sigbsq,
                betaPQL[1], betaLP[1], beta10[1], beta50[1], beta100[1],
                betaPQL[2], betaLP[2], beta10[2], beta50[2], beta100[2])

  names(estimates) <- value_labels
  return(estimates)
  
####  Might use this later for later SE / coverage analysis 
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





###############################################################
# Here we actually start running the simulation
###############################################################


# Set up some random-ish vectors to input to the generator...
# This part needs work - ask Ken - what parameters,
# and how to permute them.
n_s <- 900
p_s <- 3
b_0_s <- seq(from = -2, to = 2, by = .1)
b_1_s <- seq(from = -2, to = 2, by = .1)
b_1_s <- rev(b_1_s)
sigbsq_s <- seq(from = 1, to = 5, by = .5)

### KEY STEP - Output of the models ###
data <- as_tibble(t(mapply(check_estimates, n_s, p_s, b_0_s, b_1_s, sigbsq_s)))


# Add columns for bias of b_0 term...
# ...this could be put into the function as well.
data <- data %>% mutate(
  b_0_PQL_diff = b_0 - PQL_b_0,
  b_0_Laplace_diff = b_0 - Laplace_b_0,
  b_0_10_diff = b_0 - nAGQ_10_b_0,
  b_0_50_diff = b_0 - nAGQ_50_b_0,
  b_0_100_diff = b_0 - nAGQ_100_b_0
  )
# Add columns for bias of b_1 term...
# ...this could be put into the function as well.
data <- data %>% mutate(
  b_1_PQL_diff = b_1 - PQL_b_1,
  b_1_Laplace_diff = b_1 - Laplace_b_1,
  b_1_10_diff = b_1 - nAGQ_10_b_1,
  b_1_50_diff = b_1 - nAGQ_50_b_1,
  b_1_100_diff = b_1 - nAGQ_100_b_1
)

# Organizational structure as model list might grow
# to hopefully avoid repeated code
diff_col_names_b_0 <- c(
  "b_0_PQL_diff",
  "b_0_Laplace_diff",
  "b_0_10_diff",
  "b_0_50_diff",
  "b_0_100_diff"
  )
diff_col_names_b_1 <- c(
  "b_1_PQL_diff",
  "b_1_Laplace_diff",
  "b_1_10_diff",
  "b_1_50_diff",
  "b_1_100_diff"
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

