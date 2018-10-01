library(lme4)
library(MASS)
library(tidyverse)

# Parameters to input to the generator...
# We will explore what gives interesting results.
n <- 40
p <- 100
b_0 <- -2
b_1 <- 0.5
sigbsq <- 1

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
estimates

####  Might use this later for later SE / coverage analysis 
#  VcovPQL <- vcov(glmmresPQL, useScale = FALSE) - get covariance matrix...
# ...and use diagonal entries to get SE
#  sePQL <- sqrt(diag(VcovPQL))