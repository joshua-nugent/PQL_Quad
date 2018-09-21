library(lme4)
library(MASS)

# Vector of coefficients for the underlying model
# and constants for number of variables etc
beta <- c(-2, 1.5, 0.5, -1)
n <- 1500
p <- 3
sigbsq <- 4


# vector of p copies of each ID for 'long' data
id <- rep(1:n, each = p)

# vector of 1111...0000 for treatment arms
x1 <- as.numeric(id < (n+1)/2)

# Generate random normal values for error term, replicated p times for each ID
randint <- rep(rnorm(n, 0, sqrt(sigbsq)), each = p)

# Vector of 123...p123...p... for 'long' data
x2 <- rep(1:p, n)

# Random data from uniform distribution [0,1] for n*p cases 
x3 <- runif(p*n)

# Data generation with specified model plus noise
linpred <- beta[1] + beta[2]*x1 + beta[3]*x2 + beta[4]*x3 + randint

# Invert t
expit <- exp(linpred) / (1 + exp(linpred))

y <- runif(p*n) < expit

# ? glmmML not needed anymore?
glmmresPQL <- glmmPQL(y ~ x1 + x2 + x3, random = ~ 1 | id , family = binomial(link = "logit"))
glmmresLP <- glmer(y ~ x1 + x2 + x3 + (1|id), nAGQ = 1, family = binomial(link = "logit"))
glmmres100 <- glmer(y ~ x1 + x2 + x3 + (1|id), nAGQ = 100, family = binomial(link = "logit"))

beta
(betaPQL <- fixef(glmmresPQL))
(betaLP <- fixef(glmmresLP))
(beta100 <- fixef(glmmres100))

VcovPQL <- vcov(glmmresPQL, useScale = FALSE)
VcovLP <- vcov(glmmresLP, useScale = FALSE)
Vcov100 <- vcov(glmmres100, useScale = FALSE)

(sePQL <- sqrt(diag(VcovPQL)))
(seLP <- sqrt(diag(VcovLP)))
(se100 <- sqrt(diag(Vcov100)))


# Coverage