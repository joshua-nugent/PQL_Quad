library(lme4)
library(MASS)
library(tidyverse)
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

#sas_check <- check_estimates(n, p, b_0, b_1, sigbsq)
#sas_csv <- read.csv("data_for_sas_2.csv")
#y <- sas_csv$y
#x1 <- sas_csv$x1
#id <- sas_csv$id
#glmmresPQL <- glmmPQL(y ~ x1, random = ~ 1 | id , family = binomial(link = "logit"), niter=5)

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

  # for the SAS check we did
#  data <- cbind(id, x1, randint, linpred, expit, y)
#  write.csv(data, paste0("data_for_sas_", i,".csv"), row.names=F)
  
  glmmresPQL <- glmmPQL(y ~ x1, random = ~ 1 | id , family = binomial(link = "logit"), niter = 100)
  glmmres10 <- glmer(y ~ x1 + (1|id), nAGQ = 10, family = binomial(link = "logit"))

  betaPQL <- fixef(glmmresPQL)
  beta10 <- fixef(glmmres10)

  # get covariance matrix...
  vcovPQL <- vcov(glmmresPQL, useScale = FALSE)
  vcov10 <- vcov(glmmres10, useScale = FALSE)

  # ...and use diagonal entries to get SE
  sePQL <- sqrt(diag(vcovPQL))
  se10 <- sqrt(diag(vcov10))

  value_labels <- c("n", "p", "b_0", "b_1", "sigbsq", 
              "PQL_b_0", "nAGQ_10_b_0",
              "b_0_PQL_bias", "b_0_10_bias",
              "PQL_b_1", "nAGQ_10_b_1",
              "b_1_PQL_bias", "b_1_10_bias", 
              "b_1_PQL_SE", "b_1_10_SE"
              )
  estimates <-c(n, p, b_0, b_1, sigbsq,
                betaPQL[1], beta10[1],
                betaPQL[1] / b_0, beta10[1] / b_0,
                betaPQL[2], beta10[2],
                betaPQL[2] / b_1, beta10[2] / b_1,
                sePQL[2], se10[2]
                )

  names(estimates) <- value_labels

  return(estimates)
}


###############################################################
#
#     Here we actually start running the simulation
#
###############################################################

# Parameters for Lancet - vol 385, April 18, 2015 "Lancet 81"
# Incident Suicide attempts (Table 2, p. 1540)
# YAM vs. controls as 12 months
n <- 85  # number of clusters: 85
p <- 50 # people/observations per cluster 1987+2256=4243 people / 45+40=85 schools = 50 per cluster
b_0 <- -4.2 # intercept/baseline of control group: 1.51% = .0151 = OR of .015 = b_0 of -4.2
b_1 <- -0.8 # effect of treatment: OR of .45 = -.8
sigbsq <- 1

reps <- 100

data_attempts <- as_tibble(t(replicate(reps, check_estimates(n, p, b_0, b_1, sigbsq))))

# Parameters for Lancet - vol 385, April 18, 2015 "Lancet 81"
# Incident severe suicidal ideation: YAM vs. controls at 12 months
# Table 3, p. 1541
n <- 85  # number of clusters: 85
p <- 47 # people/observations per cluster 1991+1962=3953 people / 85 schools = 47 per cluster
b_0 <- -4.27 # intercept/baseline of control group: 1.37% = .0137 = OR of .014 = b_0 of -4.27
b_1 <- -0.69 # effect of treatment: OR of .5 = -.69
sigbsq <- 1

reps <- 100

data_ideation <- as_tibble(t(replicate(reps, check_estimates(n, p, b_0, b_1, sigbsq))))

mean(data_attempts$b_1_PQL_bias)


# # Example loops to scan over values:
# for (i in c(-2, -1.5, -1, -.5, 0, .5, 1, 1.5, 2)){
#   data <- as_tibble(t(replicate(reps, check_estimates(n, p, b_0, i, sigbsq))))
#   write.csv(data, paste0("data_for_b_1_", i,".csv"), row.names=F)
# }
# for (i in c(-1.5, -1, -.5, 0, 1, 1.5, 2)){
#   data <- as_tibble(t(replicate(reps, check_estimates(n, p, i, b_1, sigbsq))))
#   write.csv(data, paste0("data_for_b_0_", i,".csv"), row.names=F)
# }
# for (i in c(.5, 1, 1.5, 2, 3)){
#   data <- as_tibble(t(replicate(reps, check_estimates(n, p, b_0, b_1, i))))
#   write.csv(data, paste0("data_for_sigbsq_", i,".csv"), row.names=F)
# }


# Organizational structure as model list might grow
# to hopefully avoid repeated code
diff_col_names_b_0 <- c(
  "b_0_PQL_bias",
  "b_0_10_bias"
)
diff_col_names_b_1 <- c(
  "b_1_PQL_bias",
  "b_1_10_bias"
)

# read in data, make a table with bias
#b_1_bias_results <- tibble(b_1 = 0, b_1_PQL_bias = 0, b_1_Laplace_bias = 0, b_1_4_bias = 0, b_1_10_bias = 0, b_1_25_bias = 0)
b_1_bias_results <- tibble(b_1 = double(), b_1_PQL_bias = double(), b_1_Laplace_bias = double(), b_1_4_bias = double(), b_1_10_bias = double(), b_1_25_bias = double())
index <- 0
#add_row(b_1_bias_results)
for (i in c(-2, -1.5, -1, -.5, 0, .5, 1, 1.5, 2)){
  index <- index + 1
  add_row(b_1_bias_results)
  opened_file <- read.csv(paste0("data_for_b_1_", i,".csv"))
  b_1 <- i
  b_1_bias_results[index,1] <- b_1
  for (j in diff_col_names_b_1){
     m <- mean(opened_file[[j]])
    b_1_bias_results[index,j] <- m
  }
}

b_1_bias_plot <- ggplot(data = b_1_bias_results) + 
  geom_abline(slope = 0, intercept = 0) +
  geom_point(mapping = aes(x = b_1, y = b_1_PQL_bias), color = "red") +
  geom_point(mapping = aes(x = b_1, y = b_1_25_bias), color = "green")
b_1_bias_plot










##########################
############################# Below this is code junkyard
########################





 ############# integrate this into above system
# Look at bias / variance for beta1
for (i in diff_col_names_b_1){
  cat(sprintf("mean for b_1 estimate error for %s", i))
  cat(sprintf(" is: \n"))
  m <- mean(data[[i]])
  print(m)
}
for (i in diff_col_names_b_1){
  cat(sprintf("std dev for b_1 estimate error for %s", i))
  cat(sprintf(" is: \n"))
  sd <- sd(data[[i]])
  print(sd)
}
for (i in diff_col_names_b_1){
  histogram <- hist(data[[i]], main = paste("Histogram of b_1 estimate error for", i))
}

# Same as above, but for beta0 term
for (i in diff_col_names_b_0){
  cat(sprintf("mean for intercept error of %s", i))
  cat(sprintf(" is: \n"))
  m <- mean(data[[i]])
  print(m)
}
for (i in diff_col_names_b_0){
  cat(sprintf("std dev for intercept error of %s", i))
  cat(sprintf(" is: \n"))
  sd <- sd(data[[i]])
  print(sd)
}
for (i in diff_col_names_b_0){
  histogram <- hist(data[[i]], main = paste("Histogram of intercept error of", i))
}



