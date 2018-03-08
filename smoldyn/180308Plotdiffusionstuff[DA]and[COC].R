## Plotting  concentrations across the whole synapse ## 
#For file COCdiff.txt
install.packages("ggplot2")
library(ggplot2)
require(gridExtra)

### READ CSVS
molcount <- read.csv("molcount.csv");
presynapse <- read.csv("molcount_presynapse.csv");
postsynapse <- read.csv("molcount_postsynapse.csv");
cleft <- read.csv("molcount_cleft.csv");
extracellular <- read.csv("molcount_extracellular.csv");
rv1 <- read.csv("molcount_rvesicle1.csv");
rv2 <- read.csv("molcount_rvesicle2.csv");
rv3 <- read.csv("molcount_rvesicle3.csv");
rv4 <- read.csv("molcount_rvesicle4.csv");
rv5 <- read.csv("molcount_rvesicle5.csv");
rvc <- read.csv("molcount_rvesiclecentre.csv");
ru1 <- read.csv("molcount_uvesicle1.csv");
ru2 <- read.csv("molcount_uvesicle2.csv");
ru3 <- read.csv("molcount_uvesicle3.csv");
ru4 <- read.csv("molcount_uvesicle4.csv");
ru5 <- read.csv("molcount_uvesicle5.csv");
ruc <- read.csv("molcount_uvesiclecentre.csv");

time_ms <- molcount[, "time"];

#### MOLECULE LEVELS

# DAlevels 

DAoverall <- molcount[, "DA"];
DApresyn <- presynapse[, "DA"];
DAcleft <- cleft[, "DA"];
DAextracell <- extracellular[,"DA"];
DArv1 <- rv1[, "DA"];
DArv2 <- rv2[, "DA"];
DArv3 <- rv3[, "DA"];
DArv4 <- rv4[, "DA"];
DArv5 <- rv5[, "DA"];
DArvc <- rvc[, "DA"];
DAreleaesvesicles <- DArv1 + DArv2 + DArv3 + DArv4 + DArv5 + DArvc
DAuv1 <- uv1[, "DA"];
DAuv2 <- uv2[, "DA"];
DAuv3 <- uv3[, "DA"];
DAuv4 <- uv4[, "DA"];
DAuv5 <- uv5[, "DA"];
DAuvc <- uvc[, "DA"];
DAuptakevesicles <- DAuv1 + DAuv2 + DAuv3 + DAuv4 + DAuv5 + DAuvc

DAlevels <- data.frame(time_ms, DAoverall, DApresyn, DAcleft, DAextracell, 
                       DArv1, DArv2, DArv3, DArv4, DArv5, DArvc, DAreleasevesicles,
                       DAuv1, DAuv2, DAuv3, DAuv4, DAuv5, DAuvc, DAuptakevesicles);

# COC levels 

COCoverall <- molcount[, "COC"];
COCpresyn <- presynapse[, "COC"];
COCpostsyn <- postsynapse[, "COC"];
COCcleft <- cleft[, "COC"];
COCextracellular <- extracellular[, "COC"];

COClevels <- data.frame(time_ms, COCoverall, COCpresyn, COCpostsyn, COCcleft, COCextracellular)

# DAT occupancy 

DAT <- molcount[,"DAT"]
DATDA <- molcount[,"DATDA"]
DATCOC <- molcount[,"DATCOC"]

### Proportions

freeCOCcleftprop <- (COCcleft*100)/COCoverall
freeCOCextracellular <- (COCextracellular*100)/COCoverall
freeCOCprop <- (COCcleft+COCextracellular)*100/(DATCOC+COCoverall)
COCcleftprop <- (COCcleft*100)/(DATCOC+COCoverall)
COCboundprop <- (DATCOC*100)/(DATCOC+COCoverall)
COCextaprop <- (COCextracellular*100)/(DATCOC+COCoverall)

# Occupied receptors

DAToccupiedbyCOC <- (DATCOC*100)/(DAT + DATCOC + DATDA)
DAToccupiedbyDA <- (DATDA*100)/(DAT + DATDA + DATCOC)
DAToccupiedTotal <- ((DATDA+DATCOC)*100)/(DAT + DATDA + DATCOC)

Surfaceproteins <- data.frame(time_ms, DAT, DATDA, DAToccupiedbyDA, DAToccupiedbyCOC, DAToccupiedTotal)

#COC proportions

freeCOCcleftprop <- (COCcleft*100)/COCoverall
freeCOCextracellular <- (COCextracellular*100)/COCoverall
freeCOCprop <- (COCcleft+COCextracellular)*100/(DATCOC+COCoverall)
COCcleftprop <- (COCcleft*100)/(DATCOC+COCoverall)
COCboundprop <- (DATCOC*100)/(DATCOC+COCoverall)
COCextaprop <- (COCextracellular*100)/(DATCOC+COCoverall)

COCproportions <- data.frame(time_ms, freeCOCcleftprop, freeCOCextracellular, COCcleftprop,COCboundprop,COCextaprop)

## Plotting

COCproportionsfree <- ggplot(COCproportions, aes(time_ms)) +
  geom_line(aes(y = freeCOCcleftprop, colour="freeCOCcleftprop")) +
  geom_line(aes(y = freeCOCextracellular, colour="freeCOCextracellular")) +
  theme(legend.position="bottom",
      legend.text = element_text(size=9),
      legend.title = element_blank()) +
  ggtitle("free COC proportions")+
  xlab("time (ms)") +
  ylab("proportions of COC mols")

COCprops <- ggplot(COCproportions, aes(time_ms)) +
  geom_line(aes(y = freeCOCprop, colour="freeCOCprop")) +
  geom_line(aes(y = COCcleftprop, colour="COCcleftprop")) +
  geom_line(aes(y = COCboundprop, colour="COCboundprop")) +
  geom_line(aes(y = COCextaprop, colour="COCextaprop")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle(" COC proportions")+
  xlab("time (ms)") +
  ylab("proportions of COC mols")

COCnumbers <- ggplot(COClevels, aes(time_ms)) +
  geom_line(aes(y = COCoverall, colour="COCoverall")) +
  geom_line(aes(y = COCpresyn, colour="COCpresyn")) +
  geom_line(aes(y = COCpostsyn, colour="COCpostsyn")) +
  geom_line(aes(y = COCcleft, colour="COCcleft")) +
  geom_line(aes(y = COCextracellular, colour="COCextracellular")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle("COC levelsn")+
  xlab("time (ms)") +
  ylab("no. of COC mols")

DAToccupied<- ggplot(Surfaceproteins, aes(time_ms)) +
  geom_line(aes(y = DAToccupiedbyDA, colour="DA")) +
  geom_line(aes(y = DAToccupiedTotal, colour="Total")) +
  geom_line(aes(y = DAToccupiedbyCOC, colour="COC")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle("DAT occupancies")+
  xlab("time (ms)") +
  ylab("% occupied")

grid.arrange(COCnumbers,DAToccupied, COCprops, COCproportionsfree, ncol=2, nrow=2)



grid.arrange(COCnumbers,DAToccupied, COCprops, COCproportionsfree, ncol=2, nrow=2)
