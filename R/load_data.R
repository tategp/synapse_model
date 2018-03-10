#Load all data from a given experiment

#Load packages
library(dplyr)
library(ggplot2)

#Load data 
experiment <- "control"

setwd(paste0("Systems Bio/SEB/Mini Project/data/",experiment))

control <- data.frame()

for (i in 1:30) {
  fname <- paste0("molcount_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  control <- rbind(control,cbind(ID,exper,data))
}

#Filter data
control_D1DA <- control %>% 
  select(ID,time,D1DA) 
control_average_D1DA <- control_D1DA %>% 
  group_by(time) %>% 
  summarise(D1DA = mean(D1DA)) %>% 
  mutate(ID = "average")
  
ggplot(control_D1DA, aes(x = time, y = D1DA, group = ID)) +
  geom_line()
  