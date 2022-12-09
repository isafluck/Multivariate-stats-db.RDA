####### db-RDA pages(188-190 in Brocard et al.)
rm(list=ls()) 
library(vegan)
library(dplyr)

####LOAD DATA SPECIES DISTRIBUTION#####
spe_pre <- read.csv("data/PA_beetles.csv", row.names=1)
spe_pre[1:5,1:5]

spe<-spe_pre[,-1]
spe[1:5,1:5]
spe<-(spe>0)*1 #matrix of PA
colSums(spe) #all spp were recorded at some site
rowSums(spe) #all sites had spp
#no NA

#rownames siteID so the plotting gets easier
rownames(spe)<-spe_pre$siteID
spe[1:5,1:5]
dim(spe)

#####LOAD DATA ENVIRONMENT#####
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
cor(preci, temp)

####DB RDA####
#Use the "capscale" function in Vegan to run db-rda.Note that the "distance" argument turns the site by species matrix into a distance matrix. You can use any distance measure in vegan (i.e., vegdist function)

# Two ways doing db.rda --> 1) giving the raw PA matrix and specifying dist to be performed (jaccard) and function to use (vegdist)//2) perform dist function on PA before dbrda and leave other arguments empty

#####first way#####
db.rda_1 <- capscale(spe ~ elev+ long  + preci + temp , data=env, distance = "jaccard", dfun="vegdist", add=TRUE)

summary(db.rda_1)

#R2 and adjusted R2
R2 <- RsquareAdj(db.rda_1)$r.squared;R2 #21
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



#####IMPORTANCE ENVIRONMENT#####
obj<-summary(db.rda)

obj$constraints
obj$biplot#same thing as:
cor(temp, obj$constraints)
cor(preci, obj$constraints)
cor(elev, obj$constraints)
cor(long, obj$constraints)


#cummulative importance= CAP1 = 0.34 CAP2=0.61
0.61-0.34#CAP2=0.27



####COMPARE F AND Z#####
#F are site scores based on observed species distributions
#Z are site scores based on species distributions predicted by the environment.
scores(db.rda) #lc is Z, wa is F
Z<-scores(db.rda, display = c("sp","lc"), choices = 1:2)
F_values<-scores(db.rda, display = c("sp","wa"),choices = 1:2)
cor(Z$species, F_values$species) #good!

####plot Z and F scores####
x11()
par(mfrow=c(1,2))
#Fscores - observed
plot(obj$sites[,1], obj$sites[,2]) #F scores inside summary db.RDA
F_values2<-F_values$sites #extract wa using scores func
plot(db.rda,display="wa",type="n",ylab="CAP2 (27%)",xlab="CAP1 (0.34%)")
points(db.rda,display = "sites",pch=20,lwd=7)
#,col=mycolors(5)[cut(F_values2[,1],5)]

#Z scores - based env
Z2<-Z$constraints
plot(obj$constraints[,1], obj$constraints[,2], xlim=c(-1.5,2))
plot(db.rda,display="lc",type="n",ylab="CAP2 (27%)",xlab="CAP1 (0.34%)")
points(db.rda,display = "lc",pch=20,lwd=7)
#there is no "sites" option to plot the points for Z scores (because Z scores are the sites values constrained by environment/infused with environment)

####compare####
##Fscores
plot(db.rda,display="wa",type="n",ylab="CAP2 (27%)",xlab="CAP1 (0.34%)", main="Fscores")
points(db.rda,display = "sites",pch=20,lwd=7)
#Zscores
plot(db.rda,display="lc",type="n",ylab="CAP2 (27%)",xlab="CAP1 (0.34%)", main="Zscores")
points(db.rda,display = "lc",pch=20,lwd=7)

####PLOT DB.RDA####
#I'm using the first way of ferforming dbRDA because is more recent and clear aboud distances using
db.rda<-db.rda_1

#make arrows for the sites
spe.sc <- scores(db.rda, choices=1:2, scaling=1, display="sp")
site.sc <- scores(db.rda, choices=1:2, scaling=1, display="sites")


#Plot using the F-scores:####
x11()
par(mfrow=c(1,2))
plot(db.rda, scaling=1, main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores (Fscores)")
#arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=0, lty=1, col="red")
#text(site.sc, labels=rownames(spe), pos=4, cex=1)


#Plot using the Z-scores:####
plot(db.rda, scaling=1, display=c("sp","lc","cn"), main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores (Zscores)") #cn is env arrows lc display sites scores Z
#arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=0, lty=1, col="red")
#text(site.sc, labels=rownames(spe), pos=4, cex=1)



####F and Z Plots with sites names####
x11()
#Z
plot(db.rda,display="lc",type="n",ylab="CAP2 (27%)",xlab="CAP1 (0.34%)")
points(db.rda,display = "sites",pch=20,lwd=7)
#col=mycolors(5)[cut(Z$constraints[,1],5)]
text(db.rda, display="sites", col="grey30", cex=0.7, pos=3)

#F
plot(db.rda,display=c("wa"),type="n",ylab="CAP2 (27%)",xlab="CAP1 (0.34%)")
#plot(db.rda,display = "cn")
points(db.rda,display = "sites",pch=20,lwd=7)
#col=mycolors(5)[cut(Z$constraints[,1],5)]
text(db.rda, display="sites", col="grey30", cex=0.7, pos=3)



#####understanding (by hand)##### 
Yc <- as.matrix(spe)  #your site by species matrix
Xcr <- as.matrix(env)# predictors

# 2. Computation of the multivariate linear regression
# ****************************************************

# Matrix of regression coefficients (eq. 11.4)
B <- solve(t(Xcr) %*% Xcr) %*% t(Xcr) %*% Yc

# Matrix of fitted values 
Yhat <- Xcr %*% B

# Matrix of residuals
Yres <- Yc - Yhat

# Dimensions
n <- nrow(Yc)
p <- ncol(Yc)
m <- ncol(Xcr)

# 3. PCA on fitted values
# ***********************
# Covariance matrix 
S <- cov(Yhat)

# Eigenvalue decomposition-gives eigenvalues and eigenvectors (factor loadings)
eigenS <- eigen(S)

# How many canonical axes to keep/explore?
kc <- length(which(eigenS$values > 0.00000001))

# Eigenvalues of canonical axes
ev <- eigenS$values[1:kc]

# Total variance (inertia) of the response matrix
trace = sum(diag(cov(Yc)))

# Orthonormal eigenvectors (species score/FACTOR LOADINGS)
U <- eigenS$vectors[,1:kc]
row.names(U) <- colnames(Yc)

# Site scores (F-SCORES)
F <- Yc %*% U
row.names(F) <- row.names(Yc)

# Site constraints (Z-SCORES)
Z <- Yhat %*% U
row.names(Z) <- row.names(Yc)

# Canonical coefficients (regression coefficients for each explanatory variable on #each rda axis)
CC <- B %*% U
row.names(CC) <- colnames(Xcr)

# Explanatory variables
# Species-environment correlations
corXZ <- cor(Xcr,Z)

# Diagonal matrix of weights
D <- diag(sqrt(ev/trace))

# Biplot scores of explanatory variables
coordX <- corXZ %*% D    # Scaled

# Unadjusted R2
R2 <- sum(ev/trace)

# Adjusted R2
R2a <- 1-((n-1)/(n-m-1))*(1-R2)

result <- list(trace, R2, R2a, ev, CC, U, F, Z, coordX)

names(result) <- c("Total_variance", "R2", "R2adj", "Can_ev", "Can_coeff", 
                   "Species_sc1", "wa_sc1", "lc_sc1", "Biplot_sc1")

result

#compare
result
summary(db.rda)

#####ANOVA#####
#Conduct a permutation test using anova function in vegan to test the significance of the model, individual axes, and varaibles:

#Global test of the RDA result
anova(db.rda, step=1000) #the model is significant

#Tests of all canonical axes:
anova(db.rda, by="axis", step=1000) #all axis are significant

#Tests of all variables:
anova(db.rda, by="margin", step=1000)

####FOWARD MODEL SELECTION#####
step.forward <- ordiR2step(capscale(spe ~ 1, data=data.frame(env)), scope=formula(db.rda), R2scope = F, direction="forward", pstep=1000)

#the most parsimonious one:
step.forward
summary(step.forward)


######PARTIAL RDA#######
#Use the "capscale" function in Vegan to run partial db-rda

#conjuncts
db.rda <- capscale(spe ~ cbind(long+elev) + cbind(preci + temp), distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj

db.rda <- capscale(spe ~ cbind(long+elev) + Condition(preci + temp), distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj

db.rda <- capscale(spe ~ cbind(preci + temp) + Condition(long+elev), distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj


db.rda <- capscale(spe ~ long + Condition(elev + preci + temp), distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj

db.rda <- capscale(spe ~ elev + Condition(long + preci + temp), distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj

db.rda <- capscale(spe ~ temp + Condition(elev + long + preci), distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj

db.rda <- capscale(spe ~ preci + Condition(temp + elev + long), distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj

#all (regular)
db.rda <- capscale(spe ~ long+elev +preci + temp, distance = "jaccard", add=TRUE)
summary(db.rda)
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj

#####VARIANCE PARTITIONING#####
#Here we partition the variance for the model we constructed trough forward selection above: 
#first we have to create our distance matrix as the response matrix
resp<-vegdist(spe, method="jaccard")

#then we run the varpart function from vegan
#making conjuncts
spe.part <- varpart(resp,~ cbind(temp+preci),~ cbind(long+elev) ,data=env)
plot(spe.part,id.size=2,bg=3:4,Xnames=c("Env","Geo",cex=1.3))

#partition all of them
spe.part <- varpart(resp, ~ temp, ~ preci,  ~ long,  ~ elev, data=env)
spe.part
?varpart
plot(spe.part,id.size=2,bg=3:4,Xnames=c("Temp","Preci", "long", "elev",cex=1.3))

#### partial RDA x variance partitioning #####
db.rda <- capscale(spe ~ preci + Condition(temp + elev + long), distance = "jaccard", add=TRUE)
summary(db.rda) #axis 1 = 0.88
R2 <- RsquareAdj(db.rda)$r.squared;R2
R2adj <- RsquareAdj(db.rda)$adj.r.squared; R2adj #0.022

varpart(resp, ~ preci, ~ cbind(temp + elev + long), data=env) #0.023   normal=0.47

summary(db.rda) #cumulative porportion = 0.047
obj2<-scores(db.rda)
cor(preci, obj2$sites) #0.83

#The R2adj of partial RDA is similar to the individual contribution value given by varpar
