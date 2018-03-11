#Plot for DA uptake into the pre synapse 

#Load packages
library(dplyr)
library(ggplot2)

################################### Control ################################################

#Load data 
experiment <- "control"

setwd(paste0("~/Systems Bio/SEB/Mini Project/data/",experiment))

control <- data.frame()

for (i in 1:30) {
  fname <- paste0("molcount_presynapse_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  control <- rbind(control,cbind(ID,exper,data))
}

control1 <- data.frame()
for (i in 1:30) {
  fname <- paste0("molcount_releasevesicle_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  control1 <- rbind(control1,cbind(ID,exper,data))
}

#Filter data
control_DA <- as.data.frame(cbind(control$ID,control$time,control$DA,control1$DA))
colnames(control_DA) <- c("ID","time","DA","DA_release")

control_DA <- control_DA %>%
  mutate(DA = DA - DA_release) %>% 
  select(ID,time,DA)
control_average_DA <- control_DA %>% 
  group_by(time) %>% 
  summarise(DA = mean(DA)) %>% 
  mutate(ID = "average")

###########################################################################################

################################### cocaine47 ################################################

#Load data 
experiment <- "cocaine47"

setwd(paste0("~/Systems Bio/SEB/Mini Project/data/",experiment))

cocaine47 <- data.frame()
for (i in 1:30) {
  fname <- paste0("molcount_presynapse_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cocaine47 <- rbind(cocaine47,cbind(ID,exper,data))
}

cocaine471 <- data.frame()
for (i in 1:30) {
  fname <- paste0("molcount_releasevesicle_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cocaine471 <- rbind(cocaine471,cbind(ID,exper,data))
}

#Filter data
cocaine47_DA <- as.data.frame(cbind(cocaine47$ID,cocaine47$time,cocaine47$DA,cocaine471$DA))
colnames(cocaine47_DA) <- c("ID","time","DA","DA_release")

cocaine47_DA <- cocaine47_DA %>%
  mutate(DA = DA - DA_release) %>% 
  select(ID,time,DA)
cocaine47_average_DA <- cocaine47_DA %>% 
  group_by(time) %>% 
  summarise(DA = mean(DA)) %>% 
  mutate(ID = "average")

###########################################################################################

################################### cocaine77 ################################################

#Load data 
experiment <- "cocaine77"

setwd(paste0("~/Systems Bio/SEB/Mini Project/data/",experiment))

cocaine77 <- data.frame()

for (i in 1:30) {
  fname <- paste0("molcount_presynapse_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cocaine77 <- rbind(cocaine77,cbind(ID,exper,data))
}

cocaine771 <- data.frame()
for (i in 1:30) {
  fname <- paste0("molcount_releasevesicle_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cocaine771 <- rbind(cocaine771,cbind(ID,exper,data))
}

#Filter data
cocaine77_DA <- as.data.frame(cbind(cocaine77$ID,cocaine77$time,cocaine77$DA,cocaine771$DA))
colnames(cocaine77_DA) <- c("ID","time","DA","DA_release")

cocaine77_DA <- cocaine77_DA %>%
  mutate(DA = DA - DA_release) %>% 
  select(ID,time,DA)
cocaine77_average_DA <- cocaine77_DA %>% 
  group_by(time) %>% 
  summarise(DA = mean(DA)) %>% 
  mutate(ID = "average")

###########################################################################################

################################### CAD ################################################

#Load data 
experiment <- "CAD"

setwd(paste0("~/Systems Bio/SEB/Mini Project/data/",experiment))

CAD <- data.frame()

for (i in 1:30) {
  fname <- paste0("molcount_presynapse_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  CAD <- rbind(CAD,cbind(ID,exper,data))
}

CAD1 <- data.frame()
for (i in 1:30) {
  fname <- paste0("molcount_releasevesicle_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  CAD1 <- rbind(CAD1,cbind(ID,exper,data))
}

#Filter data
CAD_DA <- as.data.frame(cbind(CAD$ID,CAD$time,CAD$DA,CAD1$DA))
colnames(CAD_DA) <- c("ID","time","DA","DA_release")

CAD_DA <- CAD_DA %>%
  mutate(DA = DA - DA_release) %>% 
  select(ID,time,DA)
CAD_average_DA <- CAD_DA %>% 
  group_by(time) %>% 
  summarise(DA = mean(DA)) %>% 
  mutate(ID = "average")

###########################################################################################

################################### CCAD ################################################

#Load data 
experiment <- "CCAD"

setwd(paste0("~/Systems Bio/SEB/Mini Project/data/",experiment))

CCAD <- data.frame()

for (i in 1:30) {
  fname <- paste0("molcount_presynapse_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  CCAD <- rbind(CCAD,cbind(ID,exper,data))
}

CCAD1 <- data.frame()
for (i in 1:30) {
  fname <- paste0("molcount_releasevesicle_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  CCAD1 <- rbind(CCAD1,cbind(ID,exper,data))
}

#Filter data
CCAD_DA <- as.data.frame(cbind(CCAD$ID,CCAD$time,CCAD$DA,CCAD1$DA))
colnames(CCAD_DA) <- c("ID","time","DA","DA_release")

CCAD_DA <- CCAD_DA %>%
  mutate(DA = DA - DA_release) %>% 
  select(ID,time,DA)
CCAD_average_DA <- CCAD_DA %>% 
  group_by(time) %>% 
  summarise(DA = mean(DA)) %>% 
  mutate(ID = "average")

###########################################################################################

################################### AD ################################################

#Load data 
experiment <- "AD"

setwd(paste0("~/Systems Bio/SEB/Mini Project/data/",experiment))

AD <- data.frame()

for (i in 1:30) {
  fname <- paste0("molcount_presynapse_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  AD <- rbind(AD,cbind(ID,exper,data))
}

AD1 <- data.frame()
for (i in 1:30) {
  fname <- paste0("molcount_releasevesicle_",experiment,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  AD1 <- rbind(AD1,cbind(ID,exper,data))
}

#Filter data
AD_DA <- as.data.frame(cbind(AD$ID,AD$time,AD$DA,AD1$DA))
colnames(AD_DA) <- c("ID","time","DA","DA_release")

AD_DA <- AD_DA %>%
  mutate(DA = DA - DA_release) %>% 
  select(ID,time,DA)
AD_average_DA <- AD_DA %>% 
  group_by(time) %>% 
  summarise(DA = mean(DA)) %>% 
  mutate(ID = "average")

###########################################################################################


ggplot(data = control_DA, aes(x = time, y = DA, group = ID)) +
  geom_line(colour = "grey78") +
  geom_line(data = cocaine47_DA, aes(x = time, y = DA, group = ID), colour = "grey78") +
  geom_line(data = cocaine77_DA, aes(x = time, y = DA, group = ID), colour = "grey78") +
  geom_line(data = CAD_DA, aes(x = time, y = DA, group = ID), colour = "grey78") +
  geom_line(data = CCAD_DA, aes(x = time, y = DA, group = ID), colour = "grey78") +
  geom_line(data = AD_DA, aes(x = time, y = DA, group = ID), colour = "grey78") +
  geom_line(data = control_average_DA, aes(x = time, y = DA, colour = "control1"),size=1.5) +
  geom_line(data = cocaine47_average_DA, aes(x = time, y = DA, colour = "cocaine471"),size=1.5) +
  geom_line(data = cocaine77_average_DA, aes(x = time, y = DA, colour = "cocaine771"),size=1.5) +
  geom_line(data = CAD_average_DA, aes(x = time, y = DA, colour = "CAD1"),size=1.5) +
  geom_line(data = CCAD_average_DA, aes(x = time, y = DA, colour = "CCAD1"),size=1.5) +
  geom_line(data = AD_average_DA, aes(x = time, y = DA, colour = "AD1"),size=1.5) +
  theme_classic() + 
  scale_x_continuous(limits = c(0,250), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1000), expand = c(0, 0))
