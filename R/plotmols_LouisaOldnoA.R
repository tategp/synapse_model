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


############################# With AMPH and COCAINE #################
experiment <- "AMPHCOC"

setwd(paste0("~/SysBio - Mini-project/Smoldyn/Rev_AC_b"))

molcount_AMPHCOC <- data.frame()
cleft_AMPHCOC <- data.frame()
presynapse_AMPHCOC <- data.frame()
rvesicle_AMPHCOC <- data.frame()

for (i in 1:22) {
  fname <- paste0("molcount",2,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  molcount_AMPHCOC <- rbind(molcount_AMPHCOC,cbind(ID,exper,data))
}

for (i in 1:22) {
  fname <- paste0("molcount_cleft",2,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  cleft_AMPHCOC <- rbind(cleft_AMPHCOC,cbind(ID,exper,data))
}

for (i in 1:22) {
  fname <- paste0("molcount_presynapse",2,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  presynapse_AMPHCOC <- rbind(presynapse_AMPHCOC,cbind(ID,exper,data))
}

for (i in 1:22) {
  fname <- paste0("molcount_releasevesicle",2,"_",i,".csv")
  data <- read.table(fname,header = TRUE,sep="," )
  ID <- factor(rep(i,length(data[,1])))
  exper <- factor(rep(experiment,length(data[,1])))
  rvesicle_AMPHCOC <- rbind(rvesicle_AMPHCOC,cbind(ID,exper,data))
}

#Defining things

time_AMPHCOC <- molcount_AMPHCOC[, "time"];

#mols

DAcleft_AMPHCOC <- cleft_AMPHCOC[, "DA"];
DApresynapse_AMPHCOC <- molcount_AMPHCOC[,"DA1"]
DAvesicle_AMPHCOC <- rvesicle_AMPHCOC[,"DA"]
DAtotal_AMPHCOC <- DAcleft_AMPHCOC + DApresynapse_AMPHCOC + DAvesicle_AMPHCOC
DAtotalfree_AMPHCOC <- DAcleft_AMPHCOC + DApresynapse_AMPHCOC

AMPHpresynapse_AMPHCOC <- molcount_AMPHCOC[, "AMPH2"];
AMPHcleft_AMPHCOC <- cleft_AMPHCOC[, "AMPH"];
AMPHtotal_AMPHCOC <- AMPHpresynapse_AMPHCOC + AMPHcleft_AMPHCOC

mols_AMPHCOC <- data.frame(time_AMPHCOC, DAcleft_AMPHCOC, DApresynapse_AMPHCOC, DAvesicle_AMPHCOC,
                          DAtotal_AMPHCOC, DAtotalfree_AMPHCOC, AMPHpresynapse_AMPHCOC, AMPHcleft_AMPHCOC, 
                          AMPHtotal_AMPHCOC)

mols_AMPHCOC_average <- mols_AMPHCOC %>% 
  group_by(time_AMPHCOC) %>% 
  summarise(DAcleft_AMPHCOC = mean(DAcleft_AMPHCOC), 
            DApresynapse_AMPHCOC = mean(DApresynapse_AMPHCOC),
            DAvesicle_AMPHCOC = mean(DAvesicle_AMPHCOC), 
            DAtotal_AMPHCOC = mean(DAtotal_AMPHCOC),
            AMPHpresynapse_AMPHCOC = mean(AMPHpresynapse_AMPHCOC), 
            AMPHcleft_AMPHCOC = mean(AMPHcleft_AMPHCOC),
            AMPHtotal_AMPHCOC = mean(AMPHtotal_AMPHCOC)) %>% 
  mutate(ID = "average")

#surfaceprotiens

DAT_AMPHCOC <- molcount_AMPHCOC[, "DAT"]
DATDA_AMPHCOC <- molcount_AMPHCOC[, "DATDA"]
DATAMPH_AMPHCOC <- molcount_AMPHCOC[, "DATAMPH"] 
DATrev_AMPHCOC <- molcount_AMPHCOC[, "DATrev"]
DATrevDA_AMPHCOC <- molcount_AMPHCOC[, "DATrevDA"]
DATalltotal_AMPHCOC <- DAT_AMPHCOC + DATDA_AMPHCOC + DATAMPH_AMPHCOC + DATrev_AMPHCOC + DATrevDA_AMPHCOC
DATrevtotal_AMPHCOC <- DATrev_AMPHCOC +DATrevDA_AMPHCOC
DATtotal_AMPHCOC <- DAT_AMPHCOC + DATDA_AMPHCOC + DATAMPH_AMPHCOC
D1_AMPHCOC <- molcount_AMPHCOC[,"D1"]
D1DA_AMPHCOC <- molcount_AMPHCOC[,"D1DA"]
D1total_AMPHCOC <- D1_AMPHCOC + D1DA_AMPHCOC

surfproteins_AMPHCOC <- data.frame(time_AMPHCOC, DAT_AMPHCOC, DATDA_AMPHCOC, 
                                  DATAMPH_AMPHCOC, DATrev_AMPHCOC, DATrevDA_AMPHCOC,
                                  D1_AMPHCOC, D1DA_AMPHCOC)

surfproteins_AMPHCOC_average <- surfproteins_AMPHCOC %>% 
  group_by(time_AMPHCOC) %>% 
  summarise(DAT_AMPHCOC = mean(DAT_AMPHCOC), 
            DATDA_AMPHCOC = mean(DATDA_AMPHCOC),
            DATAMPH_AMPHCOC = mean(DATAMPH_AMPHCOC),
            DATrev_AMPHCOC = mean(DATrev_AMPHCOC),
            DATrevDA_AMPHCOC = mean(DATrevDA_AMPHCOC), 
            D1_AMPHCOC = mean(D1_AMPHCOC),
            D1DA_AMPHCOC = mean(D1DA_AMPHCOC)) %>% 
  mutate(ID = "average")


#proportions

propDAcleft_AMPHCOC <- DAcleft_AMPHCOC*100/DAtotalfree_AMPHCOC
propDApresyn_AMPHCOC <- DApresynapse_AMPHCOC*100/DAtotalfree_AMPHCOC
propAMPHcleft_AMPHCOC <- AMPHcleft_AMPHCOC*100/AMPHtotal_AMPHCOC
propAMPHpresyn_AMPHCOC <- AMPHpresynapse_AMPHCOC*100/AMPHtotal_AMPHCOC

propD1active_AMPHCOC <- D1_AMPHCOC*100/D1DA_AMPHCOC

propDAT_AMPHCOC <- DATtotal_AMPHCOC*100/DATalltotal_AMPHCOC
propDATrev_AMPHCOC <- DATrevtotal_AMPHCOC*100/DATalltotal_AMPHCOC
propDATempty_AMPHCOC <- DAT_AMPHCOC*100/DATalltotal_AMPHCOC
propDATDA_AMPHCOC <- DATDA_AMPHCOC*100/DATalltotal_AMPHCOC
propDATAMPH_AMPHCOC <- DATAMPH_AMPHCOC*100/DATalltotal_AMPHCOC
propDATrevnonactive_AMPHCOC <- DATrev_AMPHCOC*100/DATrevtotal_AMPHCOC
propDATrevDA_AMPHCOC <- DATrevDA_AMPHCOC*100/DATalltotal_AMPHCOC

proportions_AMPHCOC <- data.frame(time_AMPHCOC, propDAcleft_AMPHCOC, propDApresyn_AMPHCOC, 
                          propAMPHcleft_AMPHCOC, propAMPHpresyn_AMPHCOC, propD1active_AMPHCOC, 
                          propDAT_AMPHCOC, propDATrev_AMPHCOC, propDATempty_AMPHCOC, 
                          propDATDA_AMPHCOC, propDATAMPH_AMPHCOC, propDATrevnonactive_AMPHCOC,
                          propDATrevDA_AMPHCOC)
proportions_AMPHCOC_average <- proportions_AMPHCOC %>% 
  group_by(time_AMPHCOC) %>% 
  summarise(propDAcleft_AMPHCOC = mean(propDAcleft_AMPHCOC), 
            propDApresyn_AMPHCOC = mean(propDApresyn_AMPHCOC),
            propAMPHcleft_AMPHCOC = mean(propAMPHcleft_AMPHCOC), 
            propAMPHpresyn_AMPHCOC = mean(propAMPHpresyn_AMPHCOC),
            propD1active_AMPHCOC = mean(propD1active_AMPHCOC), 
            propDAT_AMPHCOC = mean(propDAT_AMPHCOC),
            propDATrev_AMPHCOC = mean(propDATrev_AMPHCOC),
            propDATempty_AMPHCOC = mean(propDATempty_AMPHCOC),
            propDATDA_AMPHCOC = mean(propDATDA_AMPHCOC),
            propDATAMPH_AMPHCOC = mean(propDATAMPH_AMPHCOC),
            propDATrevnonactive_AMPHCOC = mean(propDATrevnonactive_AMPHCOC),
            propDATrevDA_AMPHCOC = mean(propDATrevDA_AMPHCOC)) %>% 
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
