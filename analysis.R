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
  mutate(BH = 0.05 * seq_len(nitem) / nitem) |> 
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
  mutate(BH = 0.05 * seq_len(n_t_item) / n_t_item) |> 
  mutate(decision = p_value < BH)

Wald_results
################################################################################
# Pairwise test
################################################################################
target_items <- c(3, 5)

for(i in 1:2){
  target_item <- target_items[i]
  
  par_vec <- c()
  se_vec <- c()
  for(j in seq_along(par_est)){
    temp <- par_est[[j]][[target_item]][,2]
    par_vec <- append(par_vec, temp[1])
    se_vec <- append(se_vec, temp[2])
  }
  
  # Define group names (replace with your actual group order)
  group_names <- names(par_est)
  names(par_vec) <- group_names
  names(se_vec)  <- group_names
  
  # 1. Construct contrast matrix for all 6 pairwise comparisons (6 x 4)
  L_pair <- rbind(
    "dual vs hearing"     = c( 1, -1,  0,  0),
    "dual vs reference"   = c( 1,  0, -1,  0),
    "dual vs vision"      = c( 1,  0,  0, -1),
    "hearing vs reference"= c( 0,  1, -1,  0),
    "hearing vs vision"   = c( 0,  1,  0, -1),
    "reference vs vision" = c( 0,  0,  1, -1)
  )
  
  # 2. Compute differences and variance of differences
  diffs <- as.vector(L_pair %*% par_vec)
  var_diffs <- diag(L_pair %*% diag(se_vec^2) %*% t(L_pair))
  se_diffs  <- sqrt(var_diffs)
  
  # 3. Calculate 1-df Wald statistics and p-values
  wald_stats <- (diffs / se_diffs)^2
  p_raw      <- pchisq(wald_stats, df = 1, lower.tail = FALSE)
  
  # 4. Adjust p-values for multiple comparisons (Benjamini-Hochberg / FDR)
  p_adj_bh   <- p.adjust(p_raw, method = "BH")
  
  # 5. Combine into a clean summary data frame
  pairwise_df <- data.frame(
    Diff      = round(diffs, 3),
    SE        = round(se_diffs, 3),
    Wald_Chi2 = round(wald_stats, 3),
    p_value   = round(p_raw, 4)
  ) |> 
    arrange(p_value) |> 
    mutate(BH = 0.05 * seq_len(nrow(L_pair)) / nrow(L_pair)) |> 
    mutate(decision = p_value < BH)
  
  print(pairwise_df)
}

################################################################################
# Final Model
################################################################################
par_est <- mirt::mod2values(fit1) |> 
  filter(item %in% c("SLB038_3", "SLB038_5", "GROUP") & name %in% c("d", "MEAN_1", "COV_11")) |> 
  select(group, item, name, value)


# Reference mean and variance from the dual group
mu_ref <- par_est |>
  filter(item == "GROUP", group == "dual", name == "MEAN_1") |>
  pull(value)

var_ref <- par_est |>
  filter(item == "GROUP", group == "dual", name == "COV_11") |>
  pull(value)

# Transform group distributions
par_est <- par_est |>
  mutate(
    value = case_when(
      item == "GROUP" & name == "MEAN_1" ~
        (value - mu_ref) / sqrt(var_ref),
      
      item == "GROUP" & name == "COV_11" ~
        value / var_ref,
      
      TRUE ~ value
    )
  )

# to IRT parameter
par_est[par_est$name == "d", "value"] <- -par_est[par_est$name == "d", "value"]

par_est |> 
  filter(item == "GROUP")

par_est |> 
  filter(item != "GROUP") |> 
  arrange(item)

-coef(fit1,simplify=TRUE)$reference$items
-coef(fit1,simplify=TRUE)$reference$items[,2] |> scale() |> pnorm(lower.tail = TRUE) |> round(2)
################################################################################
# Visualization
################################################################################

library(tidyverse)

# 1. Extract group parameters (Means & SDs)
group_params <- par_est %>%
  filter(item == "GROUP") %>%
  select(group, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(sd = sqrt(COV_11))

# 2. Build curve density grid
density_data <- group_params %>%
  reframe(
    x = seq(-4, 5, length.out = 500),
    density = dnorm(x, mean = MEAN_1, sd = sd),
    .by = group
  )
density_data$group <- factor(density_data$group, levels = c("reference", "hearing", "vision", "dual"))

# 3. Filter item parameters
item_ticks_base <- par_est %>%
  filter(item != "GROUP")
item_ticks_base$item <- dplyr::recode(
  item_ticks_base$item,
  "SLB038_3" = "Q3 \n(laptop)",
  "SLB038_5" = "Q5 \n(regular cell phone)"
)
item_ticks_base$y <- 0.01
item_ticks_base$y[item_ticks_base$item == "Q3 \n(laptop)"] <- -0.01

ggplot() +
  # Curves
  geom_line(data = density_data, aes(x = x, y = density, color = group), linewidth = 1) +
  geom_ribbon(data = density_data, aes(x = x, ymin = 0, ymax = density, fill = group), alpha = 0.12) +
  # Vertical tick marks along the baseline y = 0
  geom_point(
    data = item_ticks_base,
    aes(x = value, y = y, color = group, shape = item),
    size = 5, position = position_jitter(height = 0.005, width = 0,seed = 5) # slight jitter to prevent overlap
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Digital Device Access Distributions",
    x = "Digital Device Access Level",
    y = "Density",
    color = "Group",
    fill = "Group",
    shape = "Item"
  ) +
  coord_cartesian(xlim = c(-3, 3.5), ylim = c(0, 0.6)) +
  scale_x_continuous(breaks = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

