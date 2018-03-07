## Plotting DA concentrations across the whole synapse ## 
install.packages("ggplot2")
library(ggplot2)
require(gridExtra)

molcount <- read.csv("molcount.csv");
presynapse <- read.csv("molcount_presynapse.csv");
releasevesicle <- read.csv("molcount_releasevesicle.csv");
uptakevesicle1 <- read.csv("molcount_uptakevesicle1.csv");
uptakevesicle2 <- read.csv("molcount_uptakevesicle2.csv");
#uptakevesicles <- read.csv("molcount_uptakevesicles.csv");
cleft <-read.csv("molcount_cleft.csv");

time_us <- molcount[, "time"];
time_ms <- time_us/1000

# DAlevels in 

DAoverall <- molcount[, "DA"];
DApresyn <- presynapse[, "DA"];
DAreleasevesicle <- releasevesicle[, "DA"];
DAuptakevesicle1 <- uptakevesicle1[, "DA"];
DAuptakevesicle2 <- uptakevesicle2[, "DA"];
DAuptakevesicles <- uptakevesicle1[, "DA"] + uptakevesicle2[,"DA"];
DAcleft <- (molcount[,"DA"])- (presynapse[, "DA"]);

DAlevels <- data.frame(time_ms, DAoverall, DApresyn, DAreleasevesicle, DAuptakevesicle1, DAuptakevesicle2, DAuptakevesicles, DAcleft);

#DAT and D1 occupancy with DA

DAT <- molcount[,"DAT"]
DATDA <- molcount[,"DATDA"]
DATCOC <- molcount

D1 <- molcount[,"D1"]
D1DA <- molcount[,"D1DA"]

# COC levels 

COCoverall <- molcount[, "COC"];
COCpresyn <- presynapse[, "COC"];
COCcleft <- (molcount[,"COC"])- (presynapse[, "COC"]);


DATCOC <- molcount[,"DATCOC"]

# Occupied receptors

D1occupiedbyDA <- (D1DA*100)/(D1 + D1DA)
DAToccupiedbyCOC <- (DATCOC*100)/(DAT + DATCOC + DATDA)
DAToccupiedbyDA <- (DATDA*100)/(DAT + DATDA + DATCOC)
DAToccupiedTotal <- ((DATDA+DATCOC)*100)/(DAT + DATDA + DATCOC)

Surfaceproteins <- data.frame(time_ms, DAT, DATDA, D1, D1DA, DAToccupiedbyDA, D1occupiedbyDA, DAToccupiedbyCOC, DAToccupiedTotal)

## Plotting

DAacrosssynapse <- ggplot(DAlevels, aes(time_ms)) +
  geom_line(aes(y = DAoverall, colour="DAoverall")) +
  geom_line(aes(y = DApresyn, colour="DApresyn")) +
  geom_line(aes(y = DAcleft, colour="DAcleft")) +
  theme(legend.position="bottom",
      legend.text = element_text(size=9),
      legend.title = element_blank()) +
  ggtitle("DA levels")+
  xlab("time (ms)") +
  ylab("no. of DA mols")

DAinbouton<- ggplot(DAlevels, aes(time_ms)) +
  geom_line(aes(y = DApresyn, colour="DApresyn")) +
  geom_line(aes(y = DAreleasevesicle, colour="DAreleasevesicle")) +
  geom_line(aes(y = DAuptakevesicle1, colour="DAuptakevesicle1")) +
  geom_line(aes(y = DAuptakevesicle2, colour="DAuptakevesicle2")) +
  geom_line(aes(y = DAuptakevesicles, colour="DAuptakevesicles")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle("DA levels in bouton")+
  xlab("time (ms)") +
  ylab("no. of DA mols")

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

D1occupied<- ggplot(Surfaceproteins, aes(time_ms)) +
  geom_line(aes(y = D1occupiedbyDA, colour="D1occupiedbyDA")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle("Occupied D1")+
  xlab("time (ms)") +
  ylab("% occupied")

grid.arrange(DAacrosssynapse, DAinbouton, DAToccupied, D1occupied, ncol=2, nrow=2)
