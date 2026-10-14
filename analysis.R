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

################################################################################
# NON-DIF
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

for(i in seq_len(ncol(responses))){
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
}

for(i in seq_len(ncol(responses))){
  print(anova(fit0, fits[[i]]))
}

