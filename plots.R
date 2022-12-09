#plots
library(rgdal)
library(viridis)
library(ggplot2)
mycolors<-viridis
library(vegan)

db.rda<-db.rda_1 
#load map of US
neond<-readOGR("data/NEON map/NEON_Domains.shp")
#also need objects -env- , -db.rda- from analysis script
Mean_Temperature_C<-env_pre$field_mean_annual_temperature_C


#isolate scores
smry <- summary(db.rda)
df2  <- data.frame(smry$species[,1:2])     
df3<-data.frame(smry$biplot[,1:2])
rownames(df3)<-c("Elev", "Long", "Preci", "Temp")


#######Fscores#######
#in case to scale sites
df1 <- data.frame(scores(db.rda, choices=1:2, scaling=1, display="sites"));head(df1) #more dificult read sites names
df1  <- data.frame(smry$sites[,1:2])
titulo<-"F-scores (wa)"

#######Zscores#######
df1 <- data.frame(scores(db.rda, choices=1:2, scaling=1, display="lc"));head(df1) #more difficult read sites names
df1  <- data.frame(smry$constraints[,1:2])
titulo<-"Z-scores (lc)"





#####Temperature color#####
rda.plot <- ggplot(df1, aes(x=CAP1, y=CAP2)) + 
  geom_point(aes(x=CAP1, y=CAP2, col=Mean_Temperature_C))+
  geom_text(aes(label=rownames(df1)),size=2, vjust = -1) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()+
  theme_classic()
rda.plot

####Longitude color####
Longitude<-env_pre$field_longitude
rda.plot_Z <- ggplot(df1, aes(x=CAP1, y=CAP2)) + 
  geom_point(aes(x=CAP1, y=CAP2, col=Longitude))+
  geom_text(aes(label=rownames(df1)),size=2, vjust = -1) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()+
  theme_classic()+
  scale_color_viridis_c()+
  xlim(-1.5,2)+ylim(-2,1.5)+
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df3, 
            aes(x=CAP1,y=CAP2,label=rownames(df3),
                hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
            color="red", size=3)+
  xlab("CAP1 (34%)")+
  ylab("CAP2 (27%)")+
  ggtitle(titulo)
#rda.plot


install.packages("gridExtra")
library(gridExtra)
x11()
grid.arrange(rda.plot_F, rda.plot_Z, ncol=1, nrow = 2)


#install.packages("devtools")
#devtools::install_github("gavinsimpson/ggvegan")