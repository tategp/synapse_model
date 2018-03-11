#Loading files

library(ggplot2)
require(gridExtra)
library(dplyr)

############################# Control - noAMPH #################
experiment <- "noAMPH"

setwd(paste0("~/SysBio - Mini-project/Smoldyn/Rev_noA_c"))

molcount_noAMPH <- data.frame()
cleft_noAMPH <- data.frame()
presynapse_noAMPH <- data.frame()
rvesicle_noAMPH <- data.frame()

for (i in 1:22) {
  fname <- paste0("molcount",1,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  molcount_noAMPH <- rbind(molcount_noAMPH,cbind(ID,exper,data))
}

for (i in 1:22) {
  fname <- paste0("molcount_cleft",1,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cleft_noAMPH <- rbind(cleft_noAMPH,cbind(ID,exper,data))
}

for (i in 1:22) {
  fname <- paste0("molcount_presynapse",1,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  presynapse_noAMPH <- rbind(presynapse_noAMPH,cbind(ID,exper,data))
}

for (i in 1:22) {
  fname <- paste0("molcount_releasevesicle",1,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  rvesicle_noAMPH <- rbind(rvesicle_noAMPH,cbind(ID,exper,data))
}

#Defining things

time_noAMPH <- molcount_noAMPH[, "time"];

#mols

DAcleft_noAMPH <- cleft_noAMPH[, "DA"];
DApresynapse_noAMPH <- molcount_noAMPH[,"DA1"]
DAvesicle_noAMPH <- rvesicle_noAMPH[,"DA"]
DAtotal_noAMPH <- DAcleft_noAMPH + DApresynapse_noAMPH + DAvesicle_noAMPH
DAtotalfree_noAMPH <- DAcleft_noAMPH + DApresynapse_noAMPH

AMPHpresynapse_noAMPH <- molcount_noAMPH[, "AMPH2"];
AMPHcleft_noAMPH <- cleft_noAMPH[, "AMPH"];
AMPHtotal_noAMPH <- AMPHpresynapse_noAMPH + AMPHcleft_noAMPH

mols_noAMPH <- data.frame(time_noAMPH, DAcleft_noAMPH, DApresynapse_noAMPH, DAvesicle_noAMPH,
                   DAtotal_noAMPH, DAtotalfree_noAMPH, AMPHpresynapse_noAMPH, AMPHcleft_noAMPH, 
                   AMPHtotal_noAMPH)

mols_noAMPH_average <- mols_noAMPH %>% 
  group_by(time_noAMPH) %>% 
  summarise(DAcleft_noAMPH = mean(DAcleft_noAMPH), 
            DApresynapse_noAMPH = mean(DApresynapse_noAMPH),
            DAvesicle_noAMPH = mean(DAvesicle_noAMPH), 
            DAtotal_noAMPH = mean(DAtotal_noAMPH),
            AMPHpresynapse_noAMPH = mean(AMPHpresynapse_noAMPH), 
            AMPHcleft_noAMPH = mean(AMPHcleft_noAMPH),
            AMPHtotal_noAMPH = mean(AMPHtotal_noAMPH)) %>% 
  mutate(ID = "average")

#surfaceprotiens

DAT_noAMPH <- molcount_noAMPH[, "DAT"]
DATDA_noAMPH <- molcount_noAMPH[, "DATDA"]
DATAMPH_noAMPH <- molcount_noAMPH[, "DATAMPH"] 
DATrev_noAMPH <- molcount_noAMPH[, "DATrev"]
DATrevDA_noAMPH <- molcount_noAMPH[, "DATrevDA"]
DATalltotal_noAMPH <- DAT_noAMPH + DATDA_noAMPH + DATAMPH_noAMPH + DATrev_noAMPH + DATrevDA_noAMPH
DATrevtotal_noAMPH <- DATrev_noAMPH +DATrevDA_noAMPH
DATtotal_noAMPH <- DAT_noAMPH + DATDA_noAMPH + DATAMPH_noAMPH
D1_noAMPH <- molcount_noAMPH[,"D1"]
D1DA_noAMPH <- molcount_noAMPH[,"D1DA"]
D1total_noAMPH <- D1_noAMPH + D1DA_noAMPH

surfproteins_noAMPH <- data.frame(time_noAMPH, DAT_noAMPH, DATDA_noAMPH, 
                                  DATAMPH_noAMPH, DATrev_noAMPH, DATrevDA_noAMPH,
                                  D1_noAMPH, D1DA_noAMPH)

surfproteins_noAMPH_average <- surfproteins_noAMPH %>% 
  group_by(time_noAMPH) %>% 
  summarise(DAT_noAMPH = mean(DAT_noAMPH), 
            DATDA_noAMPH = mean(DATDA_noAMPH),
            DATAMPH_noAMPH = mean(DATAMPH_noAMPH),
            DATrev_noAMPH = mean(DATrev_noAMPH),
            DATrevDA_noAMPH = mean(DATrevDA_noAMPH), 
            D1_noAMPH = mean(D1_noAMPH),
            D1DA_noAMPH = mean(D1DA_noAMPH)) %>% 
  mutate(ID = "average")

#proportions

propDAcleft_noAMPH <- DAcleft_noAMPH*100/DAtotalfree_noAMPH
propDApresyn_noAMPH <- DApresynapse_noAMPH*100/DAtotalfree_noAMPH
propAMPHcleft_noAMPH <- AMPHcleft_noAMPH*100/AMPHtotal_noAMPH
propAMPHpresyn_noAMPH <- AMPHpresynapse_noAMPH*100/AMPHtotal_noAMPH

propD1active_noAMPH <- D1_noAMPH*100/D1DA_noAMPH

propDAT_noAMPH <- DATtotal_noAMPH*100/DATalltotal_noAMPH
propDATrev_noAMPH <- DATrevtotal_noAMPH*100/DATalltotal_noAMPH
propDATempty_noAMPH <- DAT_noAMPH*100/DATalltotal_noAMPH
propDATDA_noAMPH <- DATDA_noAMPH*100/DATalltotal_noAMPH
propDATAMPH_noAMPH <- DATAMPH_noAMPH*100/DATalltotal_noAMPH
propDATrevnonactive_noAMPH <- DATrev_noAMPH*100/DATrevtotal_noAMPH
propDATrevDA_noAMPH <- DATrevDA_noAMPH*100/DATalltotal_noAMPH

proportions_noAMPH <- data.frame(time_noAMPH, 
                                 propDAcleft_noAMPH, propDApresyn_noAMPH, 
                                 propAMPHcleft_noAMPH, propAMPHpresyn_noAMPH,
                                 propD1active_noAMPH, propDAT_noAMPH, 
                                 propDATrev_noAMPH, propDATempty_noAMPH, 
                                 propDATDA_noAMPH, propDATAMPH_noAMPH,
                                 propDATrevnonactive_noAMPH, propDATrevDA_noAMPH)

proportions_noAMPH_average <- proportions_noAMPH %>% 
  group_by(time_noAMPH) %>% 
  summarise(propDAcleft_noAMPH = mean(propDAcleft_noAMPH), 
            propDApresyn_noAMPH = mean(propDApresyn_noAMPH),
            propAMPHcleft_noAMPH = mean(propAMPHcleft_noAMPH), 
            propAMPHpresyn_noAMPH = mean(propAMPHpresyn_noAMPH),
            propD1active_noAMPH = mean(propD1active_noAMPH), 
            propDAT_noAMPH = mean(propDAT_noAMPH),
            propDATrev_noAMPH = mean(propDATrev_noAMPH),
            propDATempty_noAMPH = mean(propDATempty_noAMPH),
            propDATDA_noAMPH = mean(propDATDA_noAMPH),
            propDATAMPH_noAMPH = mean(propDATAMPH_noAMPH),
            propDATrevnonactive_noAMPH = mean(propDATrevnonactive_noAMPH),
            propDATrevDA_noAMPH = mean(propDATrevDA_noAMPH)) %>% 
  mutate(ID = "average")


############################# Low AMPH #################
experiment <- "noAMPH"

setwd(paste0("~/SysBio - Mini-project/Smoldyn/Rev_noA_c"))

molcount_noAMPH <- data.frame()
cleft_noAMPH <- data.frame()

for (i in 1:22) {
  fname <- paste0("molcount",1,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  molcount_noAMPH <- rbind(molcount_noAMPH,cbind(ID,exper,data))
}

for (i in 1:22) {
  fname <- paste0("molcount_cleft",1,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cleft_noAMPH <- rbind(cleft_noAMPH,cbind(ID,exper,data))
}

#Defining things

time_noAMPH <- molcount_noAMPH[, "time"];

#mols

DAcleft_noAMPH <- cleft_noAMPH[, "DA"];
DAtotal_noAMPH <- DAcleft_noAMPH
DAtotalfree_noAMPH <- DAcleft_noAMPH + DApresynapse_noAMPH
AMPHcleft_noAMPH <- cleft_noAMPH[, "AMPH"];
AMPHtotal_noAMPH <- AMPHpresynapse_noAMPH + AMPHcleft_noAMPH

mols_noAMPH <- data.frame(time_noAMPH, DAcleft_noAMPH, DAtotal_noAMPH, 
                          DAtotalfree_noAMPH, AMPHcleft_noAMPH, 
                          AMPHtotal_noAMPH)

mols_noAMPH_average <- mols_noAMPH %>% 
  group_by(time_noAMPH) %>% 
  summarise(DAcleft_noAMPH = mean(DAcleft_noAMPH), 
            DApresynapse_noAMPH = mean(DApresynapse_noAMPH),
            DAvesicle_noAMPH = mean(DAvesicle_noAMPH), 
            DAtotal_noAMPH = mean(DAtotal_noAMPH),
            AMPHpresynapse_noAMPH = mean(AMPHpresynapse_noAMPH), 
            AMPHcleft_noAMPH = mean(AMPHcleft_noAMPH),
            AMPHtotal_noAMPH = mean(AMPHtotal_noAMPH)) %>% 
  mutate(ID = "average")

#surfaceprotiens

DAT_noAMPH <- molcount_noAMPH[, "DAT"]
DATDA_noAMPH <- molcount_noAMPH[, "DATDA"]
DATAMPH_noAMPH <- molcount_noAMPH[, "DATAMPH"] 
DATrev_noAMPH <- molcount_noAMPH[, "DATrev"]
DATrevDA_noAMPH <- molcount_noAMPH[, "DATrevDA"]
DATalltotal_noAMPH <- DAT_noAMPH + DATDA_noAMPH + DATAMPH_noAMPH + DATrev_noAMPH + DATrevDA_noAMPH
DATrevtotal_noAMPH <- DATrev_noAMPH +DATrevDA_noAMPH
DATtotal_noAMPH <- DAT_noAMPH + DATDA_noAMPH + DATAMPH_noAMPH
D1_noAMPH <- molcount_noAMPH[,"D1"]
D1DA_noAMPH <- molcount_noAMPH[,"D1DA"]
D1total_noAMPH <- D1_noAMPH + D1DA_noAMPH

surfproteins_noAMPH <- data.frame(time_noAMPH, DAT_noAMPH, DATDA_noAMPH, 
                                  DATAMPH_noAMPH, DATrev_noAMPH, DATrevDA_noAMPH,
                                  D1_noAMPH, D1DA_noAMPH)

surfproteins_noAMPH_average <- surfproteins_noAMPH %>% 
  group_by(time_noAMPH) %>% 
  summarise(DAT_noAMPH = mean(DAT_noAMPH), 
            DATDA_noAMPH = mean(DATDA_noAMPH),
            DATAMPH_noAMPH = mean(DATAMPH_noAMPH),
            DATrev_noAMPH = mean(DATrev_noAMPH),
            DATrevDA_noAMPH = mean(DATrevDA_noAMPH), 
            D1_noAMPH = mean(D1_noAMPH),
            D1DA_noAMPH = mean(D1DA_noAMPH)) %>% 
  mutate(ID = "average")

#proportions

propDAcleft_noAMPH <- DAcleft_noAMPH*100/DAtotalfree_noAMPH
propDApresyn_noAMPH <- DApresynapse_noAMPH*100/DAtotalfree_noAMPH
propAMPHcleft_noAMPH <- AMPHcleft_noAMPH*100/AMPHtotal_noAMPH
propAMPHpresyn_noAMPH <- AMPHpresynapse_noAMPH*100/AMPHtotal_noAMPH

propD1active_noAMPH <- D1_noAMPH*100/D1DA_noAMPH

propDAT_noAMPH <- DATtotal_noAMPH*100/DATalltotal_noAMPH
propDATrev_noAMPH <- DATrevtotal_noAMPH*100/DATalltotal_noAMPH
propDATempty_noAMPH <- DAT_noAMPH*100/DATalltotal_noAMPH
propDATDA_noAMPH <- DATDA_noAMPH*100/DATalltotal_noAMPH
propDATAMPH_noAMPH <- DATAMPH_noAMPH*100/DATalltotal_noAMPH
propDATrevnonactive_noAMPH <- DATrev_noAMPH*100/DATrevtotal_noAMPH
propDATrevDA_noAMPH <- DATrevDA_noAMPH*100/DATalltotal_noAMPH

proportions_noAMPH <- data.frame(time_noAMPH, 
                                 propDAcleft_noAMPH, propDApresyn_noAMPH, 
                                 propAMPHcleft_noAMPH, propAMPHpresyn_noAMPH,
                                 propD1active_noAMPH, propDAT_noAMPH, 
                                 propDATrev_noAMPH, propDATempty_noAMPH, 
                                 propDATDA_noAMPH, propDATAMPH_noAMPH,
                                 propDATrevnonactive_noAMPH, propDATrevDA_noAMPH)

proportions_noAMPH_average <- proportions_noAMPH %>% 
  group_by(time_noAMPH) %>% 
  summarise(propDAcleft_noAMPH = mean(propDAcleft_noAMPH), 
            propDApresyn_noAMPH = mean(propDApresyn_noAMPH),
            propAMPHcleft_noAMPH = mean(propAMPHcleft_noAMPH), 
            propAMPHpresyn_noAMPH = mean(propAMPHpresyn_noAMPH),
            propD1active_noAMPH = mean(propD1active_noAMPH), 
            propDAT_noAMPH = mean(propDAT_noAMPH),
            propDATrev_noAMPH = mean(propDATrev_noAMPH),
            propDATempty_noAMPH = mean(propDATempty_noAMPH),
            propDATDA_noAMPH = mean(propDATDA_noAMPH),
            propDATAMPH_noAMPH = mean(propDATAMPH_noAMPH),
            propDATrevnonactive_noAMPH = mean(propDATrevnonactive_noAMPH),
            propDATrevDA_noAMPH = mean(propDATrevDA_noAMPH)) %>% 
  mutate(ID = "average")


############################# With AMPH #################
experiment <- "withAMPH"

setwd(paste0("~/SysBio - Mini-project/Smoldyn/Rev_A_a"))

molcount_AMPH <- data.frame()
cleft_AMPH <- data.frame()
presynapse_AMPH <- data.frame()
rvesicle_AMPH <- data.frame()

for (i in 1:24) {
  fname <- paste0("molcount","_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  molcount_AMPH <- rbind(molcount_AMPH,cbind(ID,exper,data))
}

for (i in 1:24) {
  fname <- paste0("molcount_cleft","_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cleft_AMPH <- rbind(cleft_AMPH,cbind(ID,exper,data))
}

for (i in 1:24) {
  fname <- paste0("molcount_presynapse","_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  presynapse_AMPH <- rbind(presynapse_AMPH,cbind(ID,exper,data))
}

for (i in 1:24) {
  fname <- paste0("molcount_releasevesicle","_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  rvesicle_AMPH <- rbind(rvesicle_AMPH,cbind(ID,exper,data))
}

#Defining things

time_AMPH <- molcount_AMPH[, "time"];

#mols

DAcleft_AMPH <- cleft_AMPH[, "DA"];
DApresynapse_AMPH <- molcount_AMPH[,"DA1"]
DAvesicle_AMPH <- rvesicle_AMPH[,"DA"]
DAtotal_AMPH <- DAcleft_AMPH + DApresynapse_AMPH + DAvesicle_AMPH
DAtotalfree_AMPH <- DAcleft_AMPH + DApresynapse_AMPH

AMPHpresynapse_AMPH <- molcount_AMPH[, "AMPH2"];
AMPHcleft_AMPH <- cleft_AMPH[, "AMPH"];
AMPHtotal_AMPH <- AMPHpresynapse_AMPH + AMPHcleft_AMPH

mols_AMPH <- data.frame(time_AMPH, DAcleft_AMPH, DApresynapse_AMPH, DAvesicle_AMPH,
                          DAtotal_AMPH, DAtotalfree_AMPH, AMPHpresynapse_AMPH, AMPHcleft_AMPH, 
                          AMPHtotal_AMPH)

mols_AMPH_average <- mols_AMPH %>% 
  group_by(time_AMPH) %>% 
  summarise(DAcleft_AMPH = mean(DAcleft_AMPH), 
            DApresynapse_AMPH = mean(DApresynapse_AMPH),
            DAvesicle_AMPH = mean(DAvesicle_AMPH), 
            DAtotal_AMPH = mean(DAtotal_AMPH),
            AMPHpresynapse_AMPH = mean(AMPHpresynapse_AMPH), 
            AMPHcleft_AMPH = mean(AMPHcleft_AMPH),
            AMPHtotal_AMPH = mean(AMPHtotal_AMPH)) %>% 
  mutate(ID = "average")

#surfaceprotiens

DAT_AMPH <- molcount_AMPH[, "DAT"]
DATDA_AMPH <- molcount_AMPH[, "DATDA"]
DATAMPH_AMPH <- molcount_AMPH[, "DATAMPH"] 
DATrev_AMPH <- molcount_AMPH[, "DATrev"]
DATrevDA_AMPH <- molcount_AMPH[, "DATrevDA"]
DATalltotal_AMPH <- DAT_AMPH + DATDA_AMPH + DATAMPH_AMPH + DATrev_AMPH + DATrevDA_AMPH
DATrevtotal_AMPH <- DATrev_AMPH +DATrevDA_AMPH
DATtotal_AMPH <- DAT_AMPH + DATDA_AMPH + DATAMPH_AMPH
D1_AMPH <- molcount_AMPH[,"D1"]
D1DA_AMPH <- molcount_AMPH[,"D1DA"]
D1total_AMPH <- D1_AMPH + D1DA_AMPH

surfproteins_AMPH <- data.frame(time_AMPH, DAT_AMPH, DATDA_AMPH, 
                                  DATAMPH_AMPH, DATrev_AMPH, DATrevDA_AMPH,
                                  D1_AMPH, D1DA_AMPH)

surfproteins_AMPH_average <- surfproteins_AMPH %>% 
  group_by(time_AMPH) %>% 
  summarise(DAT_AMPH = mean(DAT_AMPH), 
            DATDA_AMPH = mean(DATDA_AMPH),
            DATAMPH_AMPH = mean(DATAMPH_AMPH),
            DATrev_AMPH = mean(DATrev_AMPH),
            DATrevDA_AMPH = mean(DATrevDA_AMPH), 
            D1_AMPH = mean(D1_AMPH),
            D1DA_AMPH = mean(D1DA_AMPH)) %>% 
  mutate(ID = "average")


#proportions

propDAcleft_AMPH <- DAcleft_AMPH*100/DAtotalfree_AMPH
propDApresyn_AMPH <- DApresynapse_AMPH*100/DAtotalfree_AMPH
propAMPHcleft_AMPH <- AMPHcleft_AMPH*100/AMPHtotal_AMPH
propAMPHpresyn_AMPH <- AMPHpresynapse_AMPH*100/AMPHtotal_AMPH

propD1active_AMPH <- D1_AMPH*100/D1DA_AMPH

propDAT_AMPH <- DATtotal_AMPH*100/DATalltotal_AMPH
propDATrev_AMPH <- DATrevtotal_AMPH*100/DATalltotal_AMPH
propDATempty_AMPH <- DAT_AMPH*100/DATalltotal_AMPH
propDATDA_AMPH <- DATDA_AMPH*100/DATalltotal_AMPH
propDATAMPH_AMPH <- DATAMPH_AMPH*100/DATalltotal_AMPH
propDATrevnonactive_AMPH <- DATrev_AMPH*100/DATrevtotal_AMPH
propDATrevDA_AMPH <- DATrevDA_AMPH*100/DATalltotal_AMPH

proportions_AMPH <- data.frame(time_AMPH, propDAcleft_AMPH, propDApresyn_AMPH, 
                          propAMPHcleft_AMPH, propAMPHpresyn_AMPH, propD1active_AMPH, 
                          propDAT_AMPH, propDATrev_AMPH, propDATempty_AMPH, 
                          propDATDA_AMPH, propDATAMPH_AMPH, propDATrevnonactive_AMPH,
                          propDATrevDA_AMPH)

proportions_AMPH_average <- proportions_AMPH %>% 
  group_by(time_AMPH) %>% 
  summarise(propDAcleft_AMPH = mean(propDAcleft_AMPH), 
            propDApresyn_AMPH = mean(propDApresyn_AMPH),
            propAMPHcleft_AMPH = mean(propAMPHcleft_AMPH), 
            propAMPHpresyn_AMPH = mean(propAMPHpresyn_AMPH),
            propD1active_AMPH = mean(propD1active_AMPH), 
            propDAT_AMPH = mean(propDAT_AMPH),
            propDATrev_AMPH = mean(propDATrev_AMPH),
            propDATempty_AMPH = mean(propDATempty_AMPH),
            propDATDA_AMPH = mean(propDATDA_AMPH),
            propDATAMPH_AMPH = mean(propDATAMPH_AMPH),
            propDATrevnonactive_AMPH = mean(propDATrevnonactive_AMPH),
            propDATrevDA_AMPH = mean(propDATrevDA_AMPH)) %>% 
  mutate(ID = "average")


###################################### Plotting #########

DAinpresyn <- ggplot(data = mols_noAMPH, aes(x=time_noAMPH, y = DApresynapse_noAMPH), colour = "gray78") +
  geom_line() +
  geom_line(data= mols_AMPH, aes(x=time_AMPH,y = DApresynapse_AMPH), colour =  "gray78") +
  geom_line(data= mols_noAMPH_average, aes(x=time_noAMPH,y = DApresynapse_noAMPH), colour = "orange", size = 1.2) +
  geom_line(data= mols_AMPH_average, aes(x=time_AMPH,y = DApresynapse_AMPH), colour = "red", size = 1.2) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 200), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 2500), expand = c(0,0)) +
  ggtitle("Number of DA molecules")+
  xlab("Time (ms)") +
  ylab("No. of DA in presynapse")

DAinpresyn

DAincleft <- ggplot(data= mols_noAMPH_average, aes(x=time_noAMPH,y = DAcleft_noAMPH, colour="- AMPH") +
  geom_line() +
  geom_line(data= mols_noAMPH_average, aes(x=time_noAMPH,y = DAcleft_noAMPH, colour="- AMPH"), colour = "orange", size = 1.2) +
  geom_line(data= mols_AMPH_average, aes(x=time_AMPH,y = DAcleft_AMPH, colour="+ AMPH"), colour = "red", size = 1.2) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 200), expand = c(0,0)) +
  scale_y_log10(limits = c(10, 2500), expand = c(0,0)) +
  ggtitle("Number of DA molecules")+
  xlab("Time (ms)") +
  ylab("No. of DA in Cleft")

DAincleft
