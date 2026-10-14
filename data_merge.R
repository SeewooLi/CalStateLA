library(tidyverse)

data_c <- read.csv("data/h22c_r.csv")
data_lb <- read.csv("data/h22lb_r.csv")

data_c <- data_c %>% 
  mutate(ID = as.numeric(paste(HHID, PN, sep=""))) %>% 
  select(ID, SC095, SC103) %>% 
  rename(
    vision = SC095,
    hearing = SC103) %>% 
  mutate(vision = case_when(
    vision %in% 1:3 ~ 1,
    vision %in% 4:5 ~ 0,
    .default = NA
  ),
  hearing = case_when(
    hearing %in% 1:3 ~ 1,
    hearing %in% 4:5 ~ 0,
    .default = NA
  )) %>% 
  mutate(impairment = case_when(
    vision == 1 & hearing == 1 ~ "reference",
    vision == 1 & hearing == 0 ~ "hearing",
    vision == 0 & hearing == 1 ~ "vision",
    vision == 0 & hearing == 0 ~ "dual",
    .default = NA
  ))
  
data_lb <- data_lb %>% 
  mutate(ID = as.numeric(paste(HHID, PN, sep=""))) %>% 
  select(ID, starts_with("SLB038_"))

head(data_c)
head(data_lb)

data_c <- data_c %>% 
  select(ID, impairment) %>% 
  filter(!is.na(impairment))

data_lb[data_lb==5] <- 0

################################################################################

data_merged <- left_join(data_c, data_lb, by = "ID")

responses <- data_merged %>% 
  select(starts_with("SLB038_"))

n_col <- ncol(responses)
index <- rowSums(is.na(responses)) != n_col


data_merged <- data_merged[index,]

write.csv(data_merged, "data/digital_devices.csv", row.names = FALSE)
