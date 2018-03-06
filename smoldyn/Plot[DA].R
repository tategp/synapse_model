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

time <- molcount[, "time"];

# DAlevels in 

DAoverall <- molcount[, "DA"];
DApresyn <- presynapse[, "DA"];
DAreleasevesicle <- releasevesicle[, "DA"];
DAuptakevesicle1 <- uptakevesicle1[, "DA"];
DAuptakevesicle2 <- uptakevesicle2[, "DA"];
#DAuptakevesicles <- uptakevesicles[, "DA"];
DAuptakevesicles <- uptakevesicle1[, "DA"]+uptakevesicle2[,"DA"];
DAcleft <- cleft[,"DA"];

DAlevels <- data.frame(time, DAoverall, DApresyn, DAreleasevesicle, DAuptakevesicle1, DAuptakevesicle2, DAuptakevesicles, DAcleft);

#DAT and D1 occupancy

DAT <- molcount[,"DAT"]
DATDA <- molcount[,"DATDA"]
DAToccupied <- DATDA*100/(DAT + DATDA)
D1 <- molcount[,"D1"]
D1DA <- molcount[,"D1DA"]
D1occupied <- D1DA*100/(D1 + D1DA)

Surfaceproteins <- data.frame(time, DAT, DATDA, D1, D1DA, DAToccupied, D1occupied)

## Plotting

DAacrosssynapse <- ggplot(DAlevels, aes(time)) +
  geom_line(aes(y = DAoverall, colour="DAoverall")) +
  geom_line(aes(y = DApresyn, colour="DApresyn")) +
  geom_line(aes(y = DAcleft, colour="DAcleft")) +
  theme(legend.position="bottom",
      legend.text = element_text(size=9),
      legend.title = element_blank()) +
  ggtitle("DA levels")+
  xlab("time (us)") +
  ylab("no. of DA mols")

DAinbouton<- ggplot(DAlevels, aes(time)) +
  geom_line(aes(y = DApresyn, colour="DApresyn")) +
  geom_line(aes(y = DAreleasevesicle, colour="DAreleasevesicle")) +
  geom_line(aes(y = DAuptakevesicle1, colour="DAuptakevesicle1")) +
  geom_line(aes(y = DAuptakevesicle2, colour="DAuptakevesicle2")) +
  geom_line(aes(y = DAuptakevesicles, colour="DAuptakevesicles")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle("DA levels in bouton")+
  xlab("time (us)") +
  ylab("no. of DA mols")

SPactivation<- ggplot(Surfaceproteins, aes(time)) +
  geom_line(aes(y = DAT, colour="DAT", )) +
  geom_line(aes(y = DATDA, colour="DATDA")) +
  geom_line(aes(y = D1, colour="D1")) +
  geom_line(aes(y = D1DA, colour="D1DA")) +
  scale_colour_manual("", 
                      values = c("DAT"="green", "DATDA"="darkgreen", 
                                 "D1"="blue", "D1DA"="cyan")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle("Surface protein occupancies")+
  xlab("time (us)") +
  ylab("no. of surface proteins")

SPoccupied<- ggplot(Surfaceproteins, aes(time)) +
  geom_line(aes(y = DAToccupied, colour="DAToccupied")) +
  geom_line(aes(y = D1occupied, colour="D1occupied")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  ggtitle("Surface protein occupancies")+
  xlab("time (us)") +
  ylab("% occupied")

grid.arrange(DAacrosssynapse, DAinbouton, SPactivation, SPoccupied, ncol=2, nrow=2)
