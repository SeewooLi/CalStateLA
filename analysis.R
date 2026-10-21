library(mirt)
library(tidyverse)

view_diff <- function(fit, item){
  pars <- mirt::coef(fit)
  name_groups <- names(pars)
  
  for(i in seq_along(name_groups)){
    print(name_groups[i])
    print(pars[[i]][[item]])
  }
}

################################################################################
# DATA
################################################################################
data <- read.csv("data/digital_devices.csv")

head(data)

grouping <- factor(data$impairment, levels = c("reference", "hearing", "vision", "dual"))

responses <- data[,-(1:2)]

nitem <- ncol(responses)

N <- nrow(responses)

# QUESTION:
# These next questions are about your use of modern devices. Which of the
# following devices do you own or have access to?
#
# Q1 : Desktop computer
# Q2 : iPad or other tablet
# Q3 : Laptop computer
# Q4 : Smartphone (such as iPhone, Android, or Blackberry)
# Q5 : Regular cell phone (not a smartphone)
# Q6 : E-reader (such as Kindle or Nook)
# Q7 : Wearable device (such as a Fitbit or Apple Watch)
# Q8 : Home assistant (such as Amazon Echo, Alexa, or Google Home)
# Q9 : Smart home technology or security (such as Nest, Ring, or SimpliSafe)
# Q10: Smart TV or Streaming Device (such as Roku or Amazon Firestick)



################################################################################
# NON-DIF (baseline)
################################################################################
fit0 <- mirt::multipleGroup(
  responses,
  model = 1,
  group = data$impairment, 
  itemtype = "Rasch", 
  invariance = c("free_means", 
                 "free_vars",
                 colnames(responses)
  )
)

coef(fit0)
################################################################################
# ANCHOR ITEM SELECTION (Wald 2)
################################################################################
pars <- mirt::mod2values(fit0)

mean_var_rows <- pars$item == "GROUP" 

pars[mean_var_rows,"est"] <- FALSE
pars[,"const"] <- "none"

fit1 <- mirt::multipleGroup(
  responses,
  model = 1,
  pars = pars,
  group = data$impairment, 
  itemtype = "Rasch",
  invariance = c(),
  SE = TRUE,
  SE.type = "crossprod"
)

par_est <- coef(fit1, printSE = TRUE)

Wald_results <- matrix(nrow = nitem, ncol = 4)
colnames(Wald_results) <- c("item", "X2", "df", "p_value")

for(i in seq_len(nitem)){
  cat("\n item ",i,"\n")
  view_diff(fit1, i)
  
  par_vec <- c()
  se_vec <- c()

  for(j in seq_along(par_est)){
    temp <- par_est[[j]][[i]][,2]
    par_vec <- append(par_vec, temp[1])
    se_vec <- append(se_vec, temp[2])
  }
    
  L <- rbind(
    c(1, -1,  0,  0),
    c(1,  0, -1,  0),
    c(1,  0,  0, -1)
  )
  
  diffs <- as.vector(L %*% par_vec)
  
  V <- L %*% diag(se_vec^2) %*% t(L)
  
  W <- as.numeric(t(diffs) %*% solve(V, diffs))
  df <- nrow(L)
  p_value <- pchisq(W, df = df, lower.tail = FALSE)
  
  Wald_results[i,] <- c(i, W, df, p_value)
}

(Wald_results <- as.data.frame(Wald_results))

Wald_results <- Wald_results |> 
  arrange(p_value) |> 
  mutate(BH = 0.05 / rev(seq_len(nitem))) |> 
  mutate(decision = p_value < BH)

Wald_results
################################################################################
# ANCHOR ITEM SELECTION
################################################################################

# fits <- list()
# 
# for(i in seq_len(nitem)){
#   fits[[i]] <- mirt::multipleGroup(
#     responses,
#     model = 1,
#     group = data$impairment, 
#     itemtype = "Rasch", 
#     invariance = c("free_means", 
#                    "free_vars", 
#                    colnames(responses)[-i]
#     )
#   )
#   
#   cat("\n item ",i,"\n")
#   view_diff(fits[[i]], i)
#   print(anova(fit0, fits[[i]]))
# }
# 
# 
# llik_0 <- logLik(fit0)
# for(i in seq_len(nitem)){
#   llik_a <- logLik(fits[[i]])
#   delta_llik <- llik_a - llik_0
#   
#   cat("\n item",i,": ", sqrt((delta_llik - 3)/N))
# }

################################################################################
# Wald 1
################################################################################
test_items <- Wald_results |> filter(decision) |> pull(item)
anchor_items <- setdiff(seq_len(nitem), test_items)

n_t_item <- length(test_items)

fit1 <- mirt::multipleGroup(
  responses,
  model = 1,
  group = data$impairment, 
  itemtype = "Rasch", 
  invariance = c("free_means", 
                 "free_vars",
                 colnames(responses)[anchor_items]
  ),
  SE = TRUE,
  SE.type = "crossprod"
)

par_est <- coef(fit1, printSE = TRUE)

Wald_results <- matrix(nrow = n_t_item, ncol = 4)
colnames(Wald_results) <- c("item", "X2", "df", "p_value")

for(i in seq_len(n_t_item)){
  target_item <- test_items[i]
  
  cat("\n item ",i,"\n")
  view_diff(fit1, target_item)
  
  par_vec <- c()
  se_vec <- c()
  
  for(j in seq_along(par_est)){
    temp <- par_est[[j]][[target_item]][,2]
    par_vec <- append(par_vec, temp[1])
    se_vec <- append(se_vec, temp[2])
  }
  
  L <- rbind(
    c(1, -1,  0,  0),
    c(1,  0, -1,  0),
    c(1,  0,  0, -1)
  )
  
  diffs <- as.vector(L %*% par_vec)
  
  V <- L %*% diag(se_vec^2) %*% t(L)
  
  W <- as.numeric(t(diffs) %*% solve(V, diffs))
  df <- nrow(L)
  p_value <- pchisq(W, df = df, lower.tail = FALSE)
  
  Wald_results[i,] <- c(target_item, W, df, p_value)
}

(Wald_results <- as.data.frame(Wald_results))

Wald_results <- Wald_results |> 
  arrange(p_value) |> 
  mutate(BH = 0.05 / rev(seq_len(n_t_item))) |> 
  mutate(decision = p_value < BH)

Wald_results

################################################################################
# Final Model
################################################################################
par_est <- mirt::mod2values(fit1) |> 
  filter(item %in% c("SLB038_3", "SLB038_5", "GROUP") & name %in% c("d", "MEAN_1", "COV_11")) |> 
  select(group, item, name, value)
par_est

