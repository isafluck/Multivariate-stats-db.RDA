####### db-RDA pages(188-190 in Brocard et al.)
rm(list=ls()) 
library(vegan)
library(dplyr)

setwd("C:/Users/ifluckessig/OneDrive - University of Florida/CLASSES/Fall 2022/Multivariate statistics/Lab11")

spe <- read.csv("DoubsSpe.csv", row.names=1)
head(spe)
env <- read.csv("DoubsEnv.csv", row.names=1)
head(env)

####load data from species distribution#####
setwd("C:/Users/ifluckessig/OneDrive - University of Florida/CLASSES/Fall 2022/Multivariate statistics/final project")
spe_pre <- read.csv("data from neon - beetles/PA_beetles.csv", row.names=1)
spe_pre[1:5,1:5]

spe<-spe_pre[,-1]
spe[1:5,1:5]
spe<-(spe>0)*1 #matrix of PA
colSums(spe) #all spp were recorded at some site
rowSums(spe) #all sites had spp
#no NA

####load data environment#####
env_pre<-read.csv("NEON_Field_Site_Metadata_20220412.csv", h=T) #from neon website
head(env_pre)

#filter for sites in the distribution data
index<-unique(env_pre$field_site_id)%in%unique(spe_pre$siteID)
sel_sites<-unique(env_pre$field_site_id)[index]

env_pre<-filter(env_pre,field_site_id%in%sel_sites)
head(env_pre)

#select only predictors related to hypothesis
env<-env_pre %>% 
  select(field_site_id,
         field_latitude,
         field_longitude,
         field_mean_elevation_m,
         field_mean_annual_precipitation_mm,
         field_mean_annual_temperature_C)


#geographic
lat<-env$field_latitude
long<-env$field_longitude
elev<-env$field_mean_elevation_m #I will call geographic this time because it is not related to temperature

#environment
preci<-env$field_mean_annual_precipitation_mm
temp<-env$field_mean_annual_temperature_C
#soil<-factor(env$field_soil_subgroup) - there is almost one type of soil for each site thats why I excluded

#scale
lat<-scale(lat)
long<-scale(long)
elev<-scale(elev)
preci<-scale(preci)
temp<-scale(temp)
?scale

#####db-rda didnt work#####
#Use the "capscale" function in Vegan to run db-rda.Note that the "distance" argument turns the site by species matrix into a distance matrix. You can use any distance measure in vegan (i.e., vegdist function)

db.rda <- capscale(spe ~ elev+ long  + preci + temp , data=env, distance = "jaccard", add=TRUE)

summary(db.rda)

#R2 and adjusted R2
R2 <- RsquareAdj(db.rda)$r.squared #21
R2adj <- RsquareAdj(db.rda)$adj.r.squared #13


#Plot using the F-scores:
#par(mfrow=c(1,2))
x11()
plot(db.rda, scaling=1, 
     main="Triplot db-rda F scores")
spe.sc <- scores(db.rda, choices=1:2, scaling=2, display="sp")
#arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=0, lty=1, col="red")

#Plot using the Z-scores:
plot(db.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="Triplot db-rda  Z scores")
arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=0, lty=1, col="red")

#spp in red, 
length(env_pre$field_site_id)

#####anova#####
#Conduct a permutation test using anova function in vegan to test the significance of the model, individual axes, and varaibles:

#Global test of the RDA result
anova(db.rda, step=1000)

#Tests of all canonical axes:
anova(db.rda, by="axis", step=1000)

#Tests of all variables:
anova(db.rda, by="margin", step=1000)

######partial RDA#######
#Use the "capscale" function in Vegan to run partial db-rda

#conjuncts
db.rda <- capscale(spe ~ cbind(long+elev) + cbind(preci + temp), distance = "jaccard", add=TRUE)
summary(db.rda)

#all
db.rda <- capscale(spe ~ long+elev +preci + temp, distance = "jaccard", add=TRUE)
summary(db.rda)


R2 <- RsquareAdj(db.rda)$r.squared
R2adj <- RsquareAdj(db.rda)$adj.r.squared


#Here we partition the variance for the model we constructed trough forward selection above: 
#first we have to create our distance matrix as the response matrix
resp<-vegdist(spe, method="jaccard")

#then we run the varpart function from vegan
#making conjuncts
spe.part <- varpart(resp,~ cbind(temp+preci),~ cbind(long+elev) ,data=env)
plot(spe.part,id.size=2,bg=3:4,Xnames=c("Env","Geo",cex=1.3))

#partition all of them
spe.part <- varpart(resp, ~ temp, ~ preci,  ~ long,  ~ elev, data=env)


####foward model selection#####
step.forward <- ordiR2step(capscale(spe ~ 1, data=data.frame(env2)), scope=formula(db.rda), R2scope = F, direction="forward", pstep=1000)

#the most parsimonious one:
step.forward
