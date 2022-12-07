####### db-RDA pages(188-190 in Brocard et al.)
rm(list=ls()) 
library(vegan)
library(dplyr)

####load data from species distribution#####
spe_pre <- read.csv("data/PA_beetles.csv", row.names=1)
spe_pre[1:5,1:5]

spe<-spe_pre[,-1]
spe[1:5,1:5]
spe<-(spe>0)*1 #matrix of PA
colSums(spe) #all spp were recorded at some site
rowSums(spe) #all sites had spp
#no NA

####load data environment#####
env_pre<-read.csv("data/NEON_Field_Site_Metadata_20220412.csv", h=T) #from neon website
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
env<-as.data.frame(cbind(lat,long,elev,preci,temp))
colnames(env)<- c("lat", "long", "elev", "preci", "temp")
head(env)

####db-rda####
#Use the "capscale" function in Vegan to run db-rda.Note that the "distance" argument turns the site by species matrix into a distance matrix. You can use any distance measure in vegan (i.e., vegdist function)

# Two ways doing db.rda --> 1) giving the raw PA matrix and specifying dist to be performed (jaccard) and function to use (vegdist)//2) perform dist function on PA before dbrda and leave other arguments empty

#####first way#####
db.rda_1 <- capscale(spe ~ elev+ long  + preci + temp , data=env, distance = "jaccard", dfun="vegdist", add=TRUE)

summary(db.rda_1)

#R2 and adjusted R2
R2 <- RsquareAdj(db.rda_1)$r.squared #21
R2adj <- RsquareAdj(db.rda_1)$adj.r.squared #13

#####second way#####
spe_dist<-vegdist(spe, "jaccard")
db.rda_2 <- capscale(spe_dist ~ elev+ long  + preci + temp , data=env,add=TRUE)

summary(db.rda_2)

#R2 and adjusted R2
R2 <- RsquareAdj(db.rda_2)$r.squared #21
R2adj <- RsquareAdj(db.rda_2)$adj.r.squared #13

#same results!!
#choose one
db.rda<-db.rda_1


#####check old way doing db.rda (Legendre and Gallagher 2001)#####
spe_dist<-vegdist(spe, "jaccard")#make distance of raw data PA
?cmdscale
spe_dist_pcoa<-cmdscale(spe_dist,k=5, eig=TRUE) #perform pcoa
env #predictors X

###take a look in the pcoa###
cmd<-spe_dist_pcoa
eigenvalues<-cmd$eig[1:5]
propVar<-eigenvalues/sum(eigenvalues)
cumVar<-cumsum(propVar)
PCoA_Table<-cbind(eigenvalues,propVar,cumVar)
PCoA_Table

plot(eigenvalues)
lines(lowess(eigenvalues)) #I would use 3 axis to have 70% of the variation...

#plot pcoa
x<-cmd$points[,1]
y<-cmd$points[,2]
plot(x,y,xlab= "Coordinate 1", ylab="Coordinate 2", xlim=range(x)*1.2,ylim=range(y)*1.2, type="n")
text(x,y,labels=rownames(cmd$points), cex=.9) #sites

#other way
ordiplot(scores(cmd)[,c(1,2)], type="t",cex=1, main="PCoA beetles community dissimilarity")
abline(h=0,lty=3)
abline(v=0,lty=3)

#plot spp
species<-wascores(cmd$points[,1:2],spe)
text(species,rownames(species),cex=.7, col="red") #hahah not feasible

#perform rda with the pcoa
rda_resu<-rda(spe_dist_pcoa$eig ~  elev+ long  + preci + temp , data=env)
summary(rda_resu)

R2 <- RsquareAdj(rda_resu)$r.squared #20
R2adj <- RsquareAdj(rda_resu)$adj.r.squared #12
#very similar from using db.rda function

####plot db.RDA####
#I'm using the first way of ferforming dbRDA because is more recent and clear aboud distances using
db.rda<-db.rda_1

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

#understanding
?capscale

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
