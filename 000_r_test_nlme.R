library(nlme)
library(lavaan)
library(ggplot2)
library(dplyr)


filepath <- "https://quantdev.ssri.psu.edu/sites/qdev/files/bgs_height_long.csv"
hght_long <- read.csv(file = url(filepath), header = TRUE)
head(hght_long)

summary(hght_long)

hght_long %>% group_by(id) 

# data plot
ggplot(hght_long, aes(x = age, y = hght, color = as.factor(id), group = id)) + 
  geom_point() + 
  geom_line() + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none") + 
  labs(title = "Individual Height Trajectories", y = "Height (cm)", x = "Age (months)")

# curve as estimated below - from data
b_1i = 83.08274
b_2i = 13.80617
b_3i = -3.43110
curve(b_1i + b_2i * ((x-18) / 12) + b_3i * ((x - 18) / 12) ^2, from=1, to=40, n=300, xlab="xvalue", ylab="yvalue", 
        col="blue", lwd=2)


hght.quad.nlme <- nlme(hght ~ b_1i + b_2i * ((age-18) / 12) + b_3i * ((age - 18) / 12) ^2,
                       data = hght_long,
                       fixed = b_1i + b_2i + b_3i ~ 1,
                       random =b_1i + b_2i + b_3i ~ 1,
                       groups = ~id,
                       start=c(30, 10, -3),
                       na.action=na.exclude,
                       control = nlmeControl(maxIter = 100, msVerbose = FALSE))

summary(hght.quad.nlme) 