# -------------------------------------------------------------------------------
# @project: Differential Associations of Daytime and Nighttime Heatwave Intensity 
#  with Hospitalizations for Mental and Neurological Disorders in California
# @author: Yiqun Ma
# @organization: Scripps Institution of Oceanography, UC San Diego
# @description: This script runs the case-crossover analysis to estimate the association between HWIs and hospitalizations
# @date: May 13, 2025

library(dplyr)
library(lubridate)
library(gnm)


### 1. Get stratum
data.all <- data.all %>%
  mutate(year = year(date),
         wday = wday(date),
         month = month(date),
         stratum = paste0(ZCTA,"-", wday,"-", month, "-", year))

### 2. List outcomes
outcomes <- c("Mental", "Schizophrenia", "Bipolar",
              "Depressive", "Anxiety", 
              "Conduct", "PTSD",
              "Nervous", "ADRD", "Parkinson",
              "Headache")

### 3. run for each outcome
result.list <- as.list(outcomes)
names(result.list) <- outcomes

for (j in 1:length(outcomes)){
  outcome <- outcomes[j]
  print(outcome)
  
  data <- data.all
  data[,"outcome"] <- data[,outcome]
  
  data <- data %>%
    dplyr::select(ZCTA, date, DHWI, NHWI, outcome, stratum)
  
  data <- as.data.table(data)
  data[, stratum := as.factor(stratum)]
  gc()
  
  model <- gnm(outcome ~ DHWI + NHWI,
               data = data, na.action = na.exclude,
               family = quasipoisson, eliminate = stratum)
  
  model.coef <- model$coefficients[c("DHWI","NHWI")]
  model.vcov <- vcov(model, dispersion = NULL, with.eliminate = TRUE)
  
  table <- data.frame(est = exp(model.coef),
                      SE = sqrt(c(model.vcov[1, 1], model.vcov[2, 2])))
  
  table <- table %>%
    mutate(est_lower = exp(model.coef - 1.96 * SE),
           est_upper = exp(model.coef + 1.96 * SE),
           outcome = outcome)
  
  result.list[[j]] <- table
  rm(data, model,  model.coef, model.vcov, table)
  gc()
  
}

