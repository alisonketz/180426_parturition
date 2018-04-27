###
### 4/11/2018 
### Alison Ketz
### Movement data preprocessing using adehabitatLT package
###


###
### Preliminaries
###

rm(list=ls())

library(geosphere)
library(lubridate)
library(Hmisc)
library(ggplot2)
library(adehabitatLT)
library(mvtnorm)
library(stringr)


setwd("~/Documents/Parturition/180411_parturition")

###
### load previous runs
### 

# load("anomaly_v1.Rdata")

###
### Load VIT data
###

d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.mdb',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)[10]="lowtag"
d.vit$datedropped=as.character(d.vit$datedropped)
d.vit$juliandropped=yday(mdy_hms(d.vit$datedropped))
d.vit$datefawnfound=as.character(d.vit$datefawnfound)
d.vit$julianfawnfound=yday(mdy_hms(d.vit$datefawnfound))
d.vit=d.vit[d.vit$datedropped!="",]
n.vit=length(d.vit$lowtag)


#reorder vit data by lowtag/lowtag
d.vit=d.vit[order(d.vit$lowtag),]

#extract the individual ids
individs=d.vit$lowtag

#manipulate the julian day for the doe 6865 with best guess based on movement rates

d.vit$juliandropped[9] = 147

#par(mfrow=c(1,1))
#plot(d$julian[d$lowtag==6865],d$distance[d$lowtag==6865])

###
### Load data GPS location data
###

d = matrix(NA,nr=1,nc=13)
#for loop for reading in data, using vector of lowtag's from the vit dataset
for(i in 1:n.vit){
    d.temp = read.table(paste("/home/aketz/Documents/Data/GPS_locations_ind/",d.vit$lowtag[i],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE)
    d.temp$lowtag = d.vit$lowtag[i]
    names(d.temp)=tolower(names(d.temp))
    d=rbind(d,as.matrix(d.temp))
}
d=d[-1,]
d=data.frame(d,stringsAsFactors = FALSE)

for(j in 1:dim(d)[2]){
    d[,j]=str_trim(d[,j],side="both")
}

class.indx=c(5:7,9:12)
for(j in class.indx){
    d[,j]=as.numeric(d[,j])
}

d$lowtag=as.factor(d$lowtag)

#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_local))

###
### increased fixes parturition window
###

start=yday(mdy("03/20/2017")) # May 6 beginning of parturition period
end=yday(mdy("07/07/2017")) #end of parturition period

###
### subset entire dataset to parturition window
###

d=d[d$julian>start & d$julian <= end,]

#removing last day of parturition of 5004, outlier movements
rm5004 = (dim(d[d$lowtag==5004,])[1]-15):dim(d[d$lowtag==5004,])[1]
d=d[-rm5004,]
tail(d[d$lowtag==5004,])


#or keep all data, but exclude first couple hundred observations
#for plotting all data, not just window of parturition

# d=data.frame(d %>% group_by(lowtag) %>% slice(160:(n()-300)))
# n.obs=table(d$lowtag)
# n.obs
# 
# indxing=c(1,rep(NA,n.vit-1),dim(d)[1])
# for(j in 2:n.vit){
#     indxing[j]=indxing[j-1]+n.obs[j-1]
# }
# d[indxing,]
# 
# d[d$lowtag==individs[1],]


###
### Adding Vit data to main GPS dataframe
###

d$dropped = 0

records=dim(d)[1]
records
for(i in 1:records){
    for(j in 1:n.vit){
        if(d$lowtag[i]==d.vit$lowtag[j]){
            if(d$julian[i]==d.vit$juliandropped[j])d$dropped[i]=1
        }
    }
}
sum(d$dropped)


###
### Converting date to POSIXct format
###

d$date_time_local=as.POSIXct(strptime(d$date_time_local,format="%m-%d-%Y %H:%M:%S"),tz="CST6CDT")
d$date_time_gmt=as.POSIXct(strptime(d$date_time_gmt,format="%m-%d-%Y %H:%M:%S"),tz="GMT")



#Create time lag between successive locations to censor data if needed.
time.diff <- diff(d$date_time_local)
d=d[-1,]
d$timediff <-round(as.numeric(abs(time.diff)))
rm=which(d$timediff>10)
d=d[-rm,]
names(d)[1]="devname"

###
### Dealing with missing data
###

#option 1, just remove
# d=d[!is.na(d$longitude),]

#option 2 impute
d$missing=0
for( i in 1:dim(d)[1]){
    if(is.na(d$longitude[i]))d$missing[i]=1
}

miss.per.ind=c()
for(j in 1:n.vit){
    miss.per.ind=c(miss.per.ind,sum(d$missing[d$lowtag==individs[j]]))
}
miss.per.ind
d=d[-c(1:3),]
for(i in 2:(dim(d)[1]-1)){
    if(is.na(d$longitude[i])){
        a=i-1
        while(is.na(d$longitude[a])){a=a-1}
        b=i+1
        while(is.na(d$longitude[b])){b=b+1}
        d[i,6:5] = midPoint(d[a,6:5],d[b,6:5])
    }
}

### 
### Without projection of datum into R, can use geospatial package to calculate distance and bearings
###

bearing.out=bearing(cbind(d$longitude,d$latitude))
d$bearing=bearing.out
dist.out = distHaversine(d[,6:5])
d=d[-1,]
d$distance = dist.out
d=d[-c(dim(d)[1]-1,dim(d)[1]),]#remove last 2 entries which are NA and NaN

###
### Projections! 
###

# setup coordinates
coords = cbind(d$longitude, d$latitude)
sp = SpatialPoints(coords)

# make spatial data frame
# spdf = SpatialPointsDataFrame(coords, d)
spdf = SpatialPointsDataFrame(sp, d)

# EPSG strings
latlong = "+init=epsg:4326"
proj4string(spdf) = CRS(latlong)

d.sp.proj = spTransform(spdf, CRS("+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000
                                  +y_0=-4480000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
d=data.frame(d.sp.proj)

#20 and 21 are the coordinates in UTMs x y

###
### Animal paths
### Using the adelehabitatLT

d.traj <- as.ltraj(d[,20:21], date=d$date_time_local,id=d$lowtag)
# plot(d.traj)

#Plot the trajectory for each individual - the numbers correspond to the ID in the d.traj object above
#blue is start
#red is end point

par(mfrow=c(2,2))
for(i in 1:n.vit){
    plot(d.traj[i])
}
###
### The last movements of individ1 is crazy
### removed those points up above
# plot(d.traj[1])

#converts traj object to data frame
dfdeer <- ld(d.traj)
dfdeer$id=as.character(dfdeer$id)
dfdeer=rbind(rep(NA,dim(dfdeer)[2]),dfdeer)
dfdeer=dfdeer[-dim(dfdeer)[1],]
d$rel.angle=dfdeer$rel.angle
d$dist.traj=dfdeer$dist
d$R2n=dfdeer$R2n
d$dx=dfdeer$dx
d$dy=dfdeer$dy
d$dt=dfdeer$dt
head(d)
###
### Description of parameters returned by adehabitatLT object model. 
###

# • dx, dy, dt: these parameters measured at relocation i describe the increments
# of the x and y directions and time between the relocations i and
# i + 1. Such parameters are often used in the framework of stochastic differential
# equation modelling (e.g. Brillinger et al. 2004, Wiktorsson et al.
#                     2004);
# • dist: the distance between successive relocations is often used in animal
# movement analysis (e.g. Root and Kareiva 1984, Marsh and Jones 1988);

# • abs.angle: the absolute angle αi between the x direction and the step
# built by relocations i and i + 1 is sometimes used together with the parameter
# dist to fit movement models (e.g. Marsh and Jones 1988);

# • rel.angle: the relative angle βi measures the change of direction between
# the step built by relocations i − 1 and i and the step built by relocations
# i and i + 1 (often called “turning angle”). It is often used together with
# the parameter dist to fit movement models (e.g. Root and Kareiva 1984,
#                                            Marsh and Jones 1988);

# • R2n: the squared distance between the first relocation of the trajectory
# and the current relocation is often used to test some movements models
# (e.g. the correlated random walk, see the seminal paper of Kareiva and
# Shigesada, 1983).


########################################################################################################
###
### remove the initial index of each individual 
### because it makes no sense to calculate distances 
### between the locations of the individuals
###
#######################################################################################################

remove.indx=which(is.na(d$dist.traj))
d=d[-remove.indx,]
head(d)

d$step = d$distance/d$timediff
# d$logstep = log(d$step)


########################################################################################################
###
### Calculate features
###
#######################################################################################################


###
### center and scale variables for prediction
###

scalealt=c()
scaletemp=c()
scalebear=c()
scaleangle=c()
scaleR2n=c()
scaledx=c()
scaledy=c()
scalestep=c()


#center and scale variables for features
for(j in 1:n.vit){
    scalealt = c(scalealt,scale(d$altitude[d$lowtag==individs[j]]))
    scaletemp = c(scaletemp,scale(d$temp[d$lowtag==individs[j]]))
    scalebear = c(scalebear,scale(d$bearing[d$lowtag==individs[j]]))
    scaleangle = c(scaleangle,scale(d$rel.angle[d$lowtag==individs[j]]))
    scaleR2n = c(scaleR2n,scale(d$R2n[d$lowtag==individs[j]]))
    scaledx = c(scaledx,scale(d$dx[d$lowtag==individs[j]]))
    scaledy = c(scaledy,scale(d$dy[d$lowtag==individs[j]]))
    scalestep = c(scalestep,scale(d$step[d$lowtag==individs[j]]))
}

#Test for whether I need to log transform step, and NO I don't
# Commenting out logstepscale....
# hist(logstepscale,breaks=30)
# hist(scalestep,breaks=30)
# hist(scalealt,breaks=30)

d$scalealt = scalealt
d$scaletemp = scaletemp
d$scaleangle = scaleangle
d$scalebear = scalebear
d$scaledx = scaledx
d$scaledy = scaledy
d$scalestep = scalestep
d$scaleR2n=scaleR2n

names(d)

#normal enough: altitude, temp, dx,r2n,logscalestep 
#bad: scaleangle, bearing

par(mfrow=c(3,1))
for(j in 1:n.vit){
    hist(log(d$step[d$lowtag==individs[j]]),breaks=50,main=individs[j])
    plot(density(log(d$step[d$lowtag==individs[j]])),main=individs[j])
    plot(d$julian[d$lowtag==individs[j]],log(d$step[d$lowtag==individs[j]]))
    abline(v=d.vit$juliandropped[j],col=2)
}


###
### Paturition window
###

part.window=yday("2017-05-06")

###
### Total number of possible anomalies to detect
###

possible.hits=rep(NA,n.vit)
for(j in 1:n.vit){
    possible.hits[j] = sum(d$dropped[d$lowtag==individs[j]])
}
possible.hits


##############################################################################################
###
### Expanding Anomaly Detection Training 
###
##############################################################################################
d.pre = d
d=d.pre

ph=possible.hits
names(d)
covs.train=c(13:15,30:37)
d.train = d[,covs.train]
nInd.train=length(unique(d.train[,1]))
vitdrop = d.vit$juliandropped

# source("anomalyDetectUniv.R")
# source("anomalyDetect.R")
source("training_detector.R")
source("evaluate.R")

#test 1 individual , 1 covariate
d.train = d[,covs.train]
d.train=d.train[d.train[,1]==individs[1],1:4]
epsilon=.12
fit.train=training(d.train,eps=epsilon,pw=part.window,vitdrop[1])
fit.train

hitsprob=fit.train$results_prob[fit.train$hits_indx]
outsprob=fit.train$results_prob[fit.train$outs_indx]

evaluate(alarm=fit.train$alarm,possible.hits=ph[1],nInd=1,vitdropday=vitdrop[1])
# evaluate= function(alarm,possible.hits,nInd,vitdropday){
    

#test 1 individual, multiple covariates
source("training_detector.R")
source("evaluate.R")
d.train = d[,covs.train]
d.train=d.train[d.train[,1]==individs[1],]
epsilon=rep(.12,dim(d.train)[2]-3)
fit.train=training(d.train,eps=epsilon,pw=part.window,vd = vitdrop[1])
fit.train
evaluate(alarm=fit.train$alarm,possible.hits=ph[1],nInd=1,vitdropday=d.vit$juliandropped[1])


#test multiple individuals, 1 covariate

source("training_detector.R")
source("evaluate.R")
d.train = d[,covs.train]
d.train=d.train[,c(1:3,10)]
nInd.train=length(unique(d.train[,1]))
epsilon=rep(.12,dim(d.train)[2]-3)
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop)
fit.train

evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

#test multiple individuals, multiple covariates
source("training_detector.R")
source("evaluate.R")
d.train = d[,covs.train]
epsilon=rep(.12,dim(d.train)[2]-3)
aa=Sys.time()
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop)
Sys.time()-aa
fit.train$hits_indx
evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)


###
### testing the training evaluation function
###
source("training_detector.R")
source("evaluate.R")
source("training_function.R")
d.train = d[,covs.train]

starting=Sys.time()
fit.loo=loo.train(d.train,ph=possible.hits,vdrop=vitdrop)
ending=Sys.time()
total.time=ending-starting
total.time
for(j in 1:10){
  beep(sound=4)
}


##############################################################################################
###
### Multivariate Anomaly Detection Testing
### Using all criteria
###
##############################################################################################

###precision final model
fit.prec = training(d.train,eps=fit.loo$decide[,1],pw=part.window,vd=vitdrop)# fit model
eval.prec=evaluate(alarm=fit.prec$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)

###recall final model
fit.recall = training(d.train,eps=fit.loo$decide[,2],pw=part.window,vd=vitdrop)# fit model
eval.recall=evaluate(alarm=fit.recall$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)

### F1 final model
fit.F1 = training(d.train,eps=fit.loo$decide[,3],pw=part.window,vd=vitdrop)# fit model
eval.F1=evaluate(alarm=fit.F1$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)


check.prec=rep(NA,3*n.vit)
check.recall=rep(NA,3*n.vit)
check.F1=rep(NA,3*n.vit)
check.TP=rep(NA,3*n.vit)
check.FP=rep(NA,3*n.vit)
check.FN=rep(NA,3*n.vit)


for(j in 1:(3*n.vit)){
    check.prec[j]=fit.loo$loo.eval[[j]]$out.prec
    check.recall[j]=fit.loo$loo.eval[[j]]$out.recall
    check.F1[j]=fit.loo$loo.eval[[j]]$out.F1
    check.TP[j]=fit.loo$loo.eval[[j]]$tp
    check.FP[j]=fit.loo$loo.eval[[j]]$fp
    check.FN[j]=fit.loo$loo.eval[[j]]$fn
    }
    

###recall final model
par(mfrow=c(4,1))
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly1/",individs[j],"_recall.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d$julian[d$lowtag==individs[j]],d$logscalestep[d$lowtag==individs[j]],main=individs[j])
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$scalestep,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d$julian[d$lowtag==individs[j]],d$scalealt[d$lowtag==individs[j]])
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$scalealt,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d$julian[d$lowtag==individs[j]],d$scaletemp[d$lowtag==individs[j]])
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$scaletemp,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d$julian[d$lowtag==individs[j]],d$scaledx[d$lowtag==individs[j]])
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$scaledx,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}

###
### Evaluating predictions
###

#how many days are anomalies prior to vit drop day
#col1 is total number of days with hits prior to drop+1
#col2 is the total number of days prior to vit drop
#col3 is the ratio of hits to total possible days
#col4 is the logical, whether the the algorithm hit +/- 1 day of vit drop
#col5 is the logical, whether the the algorithm hit +/- 2 days of vit dropth
#col6 is the total number of hits within 1 days
#col7 is the total number of hits within 2 days
#col8 is the total number of hits on vit drop day
#col9 is the total number of hits prior to vit drop day
hit.prec=matrix(NA,nr=n.vit,nc = 9)
for(j in 1:n.vit){
    hit.prec[j,1]=length(unique(fit.prec$alarm[[j]]$julian[fit.prec$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.prec[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.prec[j,3]=hit.prec[j,1]/hit.prec[j,2]
    hit.prec[j,4]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.prec[j,5]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.prec[j,6]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.prec[j,7]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.prec[j,8]=sum(fit.prec$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.prec[j,9]=sum(fit.prec$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
    }
hit.prec



hit.recall=matrix(NA,nr=n.vit,nc = 9)

for(j in 1:n.vit){
    hit.recall[j,1]=length(unique(fit.recall$alarm[[j]]$julian[fit.recall$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.recall[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.recall[j,3]=hit.recall[j,1]/hit.recall[j,2]
    hit.recall[j,4]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.recall[j,5]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.recall[j,6]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.recall[j,7]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.recall[j,8]=sum(fit.recall$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.recall[j,9]=sum(fit.recall$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
}
hit.recall

write.csv(hit.recall,file="summary_hit_recall_v1.csv")



hit.F1=matrix(NA,nr=n.vit,nc = 9)
for(j in 1:n.vit){
  hit.F1[j,1]=length(unique(fit.F1$alarm[[j]]$julian[fit.F1$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
  hit.F1[j,2]=(d.vit$juliandropped[j]-1)-126+1
  hit.F1[j,3]=hit.F1[j,1]/hit.F1[j,2]
  hit.F1[j,4]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
  hit.F1[j,5]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
  hit.F1[j,6]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
  hit.F1[j,7]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
  hit.F1[j,8]=sum(fit.F1$alarm[[j]]$julian ==d.vit$juliandropped[j])
  hit.F1[j,9]=sum(fit.F1$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
}
hit.F1

write.csv(hit.F1,file="summary_hit_F1_v1.csv")

### output for recall

Metric=c("total number of days with hits prior to day before vit drop",
         "total number of days prior to vit drop",
         "ratio of hits to total possible days",
         "logical- whether the the algorithm hit +/- 1 day of vit drop",
         "logical- whether the the algorithm hit +/- 2 days of vit dropth",
         "total number of hits within 1 days",
         "total number of hits within 2 days",
         "total number of hits on vit drop day",
         "total number of hits prior to vit drop day")

hit.recall[,3]=round(hit.recall[,3],2)

out.recall=as.data.frame(cbind(Metric,t(hit.recall)),stringsAsFactors=FALSE)
names(out.recall)=c("Metric",as.character(individs))
library(xtable)
xtable(out.recall,include.rownames=FALSE)



###
### population level zscores for step and alt
###

d$zstep = scale(d$step)
d$zalt = scale(d$altitude)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pstepz = ggplot(d,aes(d$zstep,color=lowtag))+geom_density(size=1.01,alpha=.5)+xlim(-3,5)+theme_bw()+
    scale_color_manual(values=c(rep("grey70",8),cbPalette[2],cbPalette[3],cbPalette[4],"grey70"))+ylab("Density")+xlab("Step rate")
suppressWarnings(pstepz)

paltz = ggplot(d,aes(d$zalt,color=lowtag))+geom_density(size=1.01,alpha=.5)+xlim(-3,5)+theme_bw()+
    scale_color_manual(values=c(rep("grey70",8),cbPalette[2],cbPalette[3],cbPalette[4],"grey70"))
paltz


save.image("anomaly_v1.Rdata")
