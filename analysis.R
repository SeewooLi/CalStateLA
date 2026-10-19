library(mirt)

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
# ANCHOR ITEM SELECTION
################################################################################
fits <- list()

for(i in seq_len(nitem)){
  fits[[i]] <- mirt::multipleGroup(
    responses,
    model = 1,
    group = data$impairment, 
    itemtype = "Rasch", 
    invariance = c("free_means", 
                   "free_vars", 
                   colnames(responses)[-i]
    )
  )
  
  cat("\n item ",i,"\n")
  view_diff(fits[[i]], i)
  print(anova(fit0, fits[[i]]))
}


llik_0 <- logLik(fit0)
for(i in seq_len(nitem)){
  llik_a <- logLik(fits[[i]])
  delta_llik <- llik_a - llik_0
  
  cat("\n item",i,": ", sqrt((delta_llik - 3)/N))
}

################################################################################
# Second-stage Filtering
################################################################################
anchor_items <- c(2,6,8,9,10)
test_items <- setdiff(seq_len(nitem), anchor_items)

fit0 <- mirt::multipleGroup(
  responses,
  model = 1,
  group = data$impairment, 
  itemtype = "Rasch", 
  invariance = c("free_means", 
                 "free_vars",
                 colnames(responses)[anchor_items]
  )
)

fits <- list()

for(i in seq_along(test_items)){
  fits[[i]] <- mirt::multipleGroup(
    responses,
    model = 1,
    group = data$impairment, 
    itemtype = "Rasch", 
    invariance = c("free_means", 
                   "free_vars", 
                   colnames(responses)[c(anchor_items, test_items[i])]
    )
  )
  
  cat("\n item ",test_items[i],"\n")
  view_diff(fits[[i]], i)
  print(anova(fits[[i]], fit0))
}


llik_0 <- logLik(fit0)
for(i in seq_along(test_items)){
  llik_a <- logLik(fits[[i]])
  delta_llik <- llik_0 - llik_a
  
  cat("\n item",test_items[i],": ", sqrt((delta_llik - 3)/N))
}

################################################################################
# Third-stage Filtering
################################################################################
anchor_items <- c(2,6,7,8,9,10)
test_items <- setdiff(seq_len(nitem), anchor_items)

fit0 <- mirt::multipleGroup(
  responses,
  model = 1,
  group = data$impairment, 
  itemtype = "Rasch", 
  invariance = c("free_means", 
                 "free_vars",
                 colnames(responses)[anchor_items]
  )
)

fits <- list()

for(i in seq_along(test_items)){
  fits[[i]] <- mirt::multipleGroup(
    responses,
    model = 1,
    group = data$impairment, 
    itemtype = "Rasch", 
    invariance = c("free_means", 
                   "free_vars", 
                   colnames(responses)[c(anchor_items, test_items[i])]
    )
  )
  
  cat("\n item ",test_items[i],"\n")
  view_diff(fits[[i]], test_items[i])
  print(anova(fits[[i]], fit0))
}


llik_0 <- logLik(fit0)
for(i in seq_along(test_items)){
  llik_a <- logLik(fits[[i]])
  delta_llik <- llik_0 - llik_a
  
  cat("\n item",test_items[i],": ", sqrt((delta_llik - 3)/N))
}

################################################################################
# Final Model
################################################################################
for(i in seq_along(test_items)){
  cat("\n item ",test_items[i],"\n")
  view_diff(fit0, test_items[i])
}



