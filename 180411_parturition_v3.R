##################################################################################################
###
### 4/11/2018 
### Alison Ketz
### Anomaly detection training using movement data with data processing
###
##################################################################################################

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
library(beepr)
library(stringr)

setwd("~/Documents/Parturition/180426_parturition")

###
### Load previous runs
###

# load("anomaly_v3.Rdata")

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

# d.vit=d.vit[-9,]
# individs=d.vit$lowtag
# n.vit=length(d.vit$lowtag)

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
### Projections
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

#Plot the trajectory for each individual - the numbers correspond to the ID in the d.traj object above
#blue is start
#red is end point
# par(mfrow=c(2,2))
# for(i in 1:n.vit){
#     plot(d.traj[i])
# }

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


remove.indx=which(is.na(d$dist.traj))
d=d[-remove.indx,]

d$step = d$distance/d$timediff
# d$logstep = log(d$step)

##############################################################################################
###
### creating a two day based time step data frame, and running anomaly detection 
###
##############################################################################################
# features to calculate within day means and standard deviations, and assign to daily
#lowtag, julian, dropped,step,altittude, temp, bearing, rel.angle, R2n, dx,dy
#lowtag, julian, dropped,
# mean and sd of : step,altitude, temp, bearing, rel.angle, R2n, dx,dy
#matrix dim should  = 8*2+3 = 19

d$lowtag = as.character(d$lowtag)
d.dos.temp=matrix(NA,nr = 1,nc=19)#this will be too big, will trim after
for(j in 1:n.vit){
    d.temp=d[d$lowtag == individs[j],]
    julian.temp = unique(d.temp$julian)
    temp.mat = matrix(NA,nr = length(julian.temp),nc = 19)
    for(i in 1:length(julian.temp)){
        temp.mat[i,1] = d.temp$lowtag[1]
        temp.mat[i,2] = julian.temp[i]
        temp.mat[i,3] = ifelse(sum(d.temp$dropped[d.temp$julian == julian.temp[i]])>1,1,0)
        temp.mat[i,4] = mean(d.temp$step[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.step
        temp.mat[i,5] = sd(d.temp$step[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.step
        temp.mat[i,6] = mean(d.temp$altitude[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.altitude
        temp.mat[i,7] = sd(d.temp$altitude[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.altitude
        temp.mat[i,8] = mean(d.temp$temp[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.temp
        temp.mat[i,9] = sd(d.temp$temp[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.temp
        temp.mat[i,10] = mean(d.temp$bearing[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.bearing
        temp.mat[i,11] = sd(d.temp$bearing[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.bearing
        temp.mat[i,12] = mean(d.temp$rel.angle[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.rel.angle
        temp.mat[i,13] = sd(d.temp$rel.angle[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.rel.angle
        temp.mat[i,14] = mean(d.temp$R2n[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.R2n
        temp.mat[i,15] = sd(d.temp$R2n[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.R2n
        temp.mat[i,16] = mean(d.temp$dx[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dx
        temp.mat[i,17] = sd(d.temp$dx[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dx
        temp.mat[i,18] = mean(d.temp$dy[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dy
        temp.mat[i,19] = sd(d.temp$dy[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dy
    }
    d.dos.temp = rbind(d.dos.temp,temp.mat)
}
d.dos.temp=d.dos.temp[-1,]
d.dos = data.frame(d.dos.temp,stringsAsFactors = FALSE)
names(d.dos)=c("lowtag","julian","dropped","mu.step","sig.step","mu.altitude","sig.altitude","mu.temp","sig.temp",
                "mu.bearing","sig.bearing","mu.rel.angle","sig.rel.angle","mu.R2n","sig.R2n","mu.dx","sig.dx","mu.dy","sig.dy")
for(i in 3:19){
    d.dos[,i]=as.numeric(d.dos[,i])
}

###
### center and scale variables for prediction
###

for(j in 1:n.vit){
    d.dos[d.dos$lowtag == individs[j],4:19]=apply(d.dos[d.dos$lowtag == individs[j],4:19],2,scale)
}
d.dos$julian = as.integer(d.dos$julian)

head(d.dos)

###
### Paturition window
###

part.window=yday("2017-05-06")

###
### Total number of possible anomalies to detect
###

possible.hits=rep(NA,n.vit)
for(j in 1:n.vit){
    possible.hits[j] = sum(d.dos$dropped[d.dos$lowtag==individs[j]])
}

##############################################################################################
###
### Debugging all cases of Anomaly Detection Training 
###
##############################################################################################

d.train=d.dos
ph=possible.hits
nInd.train=length(unique(d.train[,1]))
vitdrop = d.vit$juliandropped

source("training_detector.R")
source("training_test_detector.R")
source("evaluate.R")

#test 1 individual , 1 covariate
source("training_detector.R")
source("evaluate.R")
d.train=d.train[d.train[,1]==individs[1],1:4]
epsilon=.15
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop[1])
fit.train$alarm[[1]]
evaluate(alarm=fit.train$alarm,possible.hits=ph[1],nInd=1,vitdropday=vitdrop[1])
fit.train.test=training_test(d.train,eps=epsilon,pw=part.window,vd=vitdrop[1])

# nCovs=1
# fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs[1],d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model
# evaluate(alarm=fit.train2$alarm,possible.hits=ph[1],nInd=1,vitdropday=vitdrop[1])
# fit.train$ind.mean
# fit.train$ind.sd
# fit.train2$ind.mean
# fit.train2$ind.sd
# fit.train$detect.quant
# fit.train2$detect.quant
# fit.train$threshold.density
# fit.train2$threshold.density
# fit.train$detect.density
# fit.train2$detect.density
# fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]

#test 1 individual, multiple covariates
# source("anomalyDetect.R")
source("training_detector.R")
source("training_test_detector.R")
source("evaluate.R")

d.train = d.dos
d.train=d.train[d.train[,1]==individs[1],]
epsilon=rep(.12,dim(d.train)[2]-3)
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop[1])
fit.train
evaluate(alarm=fit.train$alarm,possible.hits=ph[1],nInd=1,vitdropday=vitdrop[1])
training_test(d.train,eps=epsilon,pw=part.window,vd=vitdrop[1])

# nCovs=16
# fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs[1],d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model
# evaluate(alarm=fit.train2$alarm,possible.hits=ph[1],nInd=1,vitdropday=vitdrop[1])
# fit.train$ind.mean==fit.train2$ind.mean
# fit.train$ind.sd==fit.train2$ind.sd
# fit.train$detect.quant
# fit.train2$detect.quant
# fit.train$threshold.density
# fit.train2$threshold.density
# fit.train$detect.density
# fit.train2$detect.density
# fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]
# fit.train$alarm[[1]][,2]
# fit.train2$alarm[[1]][,2]

#test multiple individuals, 1 covariate
source("training_detector.R")
source("training_test_detector.R")
source("evaluate.R")

d.train = d.dos
d.train=d.train[,1:4]
nInd.train=length(unique(d.train[,1]))
epsilon=rep(.12,dim(d.train)[2]-3)
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop)
fit.train$alarm
evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=12,vitdropday=vitdrop)
training_test(d.train,eps=epsilon,pw=part.window,vd=vitdrop)

# nCovs=1
# fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs,d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model
# evaluate(alarm=fit.train2$alarm,possible.hits=ph,nInd=12,vitdropday=vitdrop)
# fit.train$ind.mean==fit.train2$ind.mean
# fit.train$ind.mean
# fit.train2$ind.mean
# fit.train$ind.sd==fit.train2$ind.sd
# fit.train$detect.quant
# fit.train2$detect.quant
# fit.train$threshold.density
# fit.train2$threshold.density
# fit.train$detect.density
# fit.train2$detect.density
# fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]
# fit.train$alarm[[1]][,2]
# fit.train2$alarm[[1]][,2]


#test multiple individuals, multiple covariates
source("training_detector.R")
source("training_test_detector.R")
source("evaluate.R")

d.train = d.dos
epsilon=rep(.15,dim(d.train)[2]-3)
aa=Sys.time()
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop)
Sys.time()-aa
fit.train$alarm[[1]]
evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=12,vitdropday=vitdrop)
training_test(d.train,eps=epsilon,pw=part.window,vd=vitdrop)

# nCovs=16
# fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs,d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model
# evaluate(alarm=fit.train2$alarm,possible.hits=ph,nInd=12,vitdropday=vitdrop)
# fit.train$ind.mean==fit.train2$ind.mean
# fit.train$ind.mean
# fit.train2$ind.mean
# fit.train$ind.sd==fit.train2$ind.sd
# fit.train$detect.quant
# fit.train2$detect.quant
# fit.train$threshold.density
# fit.train2$threshold.density
# fit.train$detect.density
# fit.train2$detect.density
# fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]
# fit.train$alarm[[1]][,2]
# fit.train2$alarm[[1]][,2]

##############################################################################################
##############################################################################################
###
### Anomaly Detection Training 
###
##############################################################################################
##############################################################################################

source("training_function_v1.R")
source("training_function_v2.R")
source("training_function_v3.R")

source("training_detector.R")
source("training_test_detector.R")
source("evaluate.R")
# source("anomalyDetect.R")

# starting=Sys.time()
fit.loo.v1=loo.train.v1(d.train=d.dos,part.window=126,ph=possible.hits,vdrop=vitdrop)
# ending=Sys.time()
# total.time = ending-starting
# total.time

# starting=Sys.time()
fit.loo.v2=loo.train.v2(d.train=d.dos,part.window=126,ph=possible.hits,vdrop=vitdrop)
# ending=Sys.time()
# total.time = ending-starting
# total.time

fit.loo=fit.loo.v1

feature.indx=rep(1:16,2)
version=c(rep("v1",16),rep("v2",16))
prec=c(fit.loo.v1$decide.star[,1],fit.loo.v2$decide[,1])
recall=c(fit.loo.v1$decide.star[,2],fit.loo.v2$decide[,2])
F1=c(fit.loo.v1$decide.star[,3],fit.loo.v2$decide[,3])

test.df =data.frame(feature.indx,version,prec,recall,F1)
test.df
pdf("Test_training_versions_m6.pdf")
ggplot(test.df,aes(x=feature.indx))+geom_point(aes(y=prec,color=version),alpha=.5,size=3)+scale_color_manual(values=c("#E69F00", "#56B4E9"))+theme_bw()
dev.off()

for(j in 1:10){
  beep(sound=5)
}


#fitting with different data for determining the quantiles for the density threshold calculation
# fit.loo.v1.m1 = fit.loo.v1#<(vd-2),test <(vd-2)
# fit.loo.v1.m2 = fit.loo.v1#=1:nT,test <(vd-2)
# fit.loo.v1.m3 = fit.loo.v1#<pw,test <(vd-2)
# fit.loo.v1.m4 = fit.loo.v1#>=pw,test <(vd-2)
# fit.loo.v1.m5 = fit.loo.v1#<(vd-2), test <pw
fit.loo.v1.m6 = fit.loo.v1#=1:nT,test <pw
# fit.loo.v1.m7 = fit.loo.v1#<pw,test <pw
# fit.loo.v1.m8 = fit.loo.v1#>=pw,test <pw

###
### v1.m1,v1.m2,v1.m3,v1.m4
### v1.m5,v1.m6,v1.m7,v2.m8

#fitting with different data for determining the quantiles for the density threshold calculation
# fit.loo.v2.m1 = fit.loo.v2#<(vd-2),test <(vd-2)
# fit.loo.v2.m2 = fit.loo.v2#=1:nT,test <(vd-2)
# fit.loo.v2.m3 = fit.loo.v2#<pw,test <(vd-2)
# fit.loo.v2.m4 = fit.loo.v2#>=pw,test <(vd-2)
# fit.loo.v2.m5 = fit.loo.v2#<(vd-2), test <pw
fit.loo.v2.m6 = fit.loo.v2#=1:nT,test <pw
# fit.loo.v2.m7 = fit.loo.v2#<pw,test <pw
# fit.loo.v2.m8 = fit.loo.v2#>=pw,test <pw

###
### v2.m1,v2.m2,v2.m3,v2.m4,
### v2.m5,v2.m6,v2.m7,v2.m8
###


# loo_train = function(d.train,part.window=126,ph,vdrop){
# ending=Sys.time()
# total.time=ending-starting
# total.time


attributes(fit.loo.v1.m1)


fit.loo$decide
fit.loo$epsilon.star

fit.loo$train.cov[,,2]
fit.loo$compare
fit.loo$loo.eval
fit.loo$train.cov[,,2]
fit.loo$loo.fittest


##############################################################################################
###
### Evaluate model training
### 
###
##############################################################################################


recall.v1 = c(mean(fit.loo.v1.m1$eval[,2,2]),
  mean(fit.loo.v1.m2$eval[,2,2]),
  mean(fit.loo.v1.m3$eval[,2,2]),
  mean(fit.loo.v1.m4$eval[,2,2]),
  mean(fit.loo.v1.m5$eval[,2,2]),
  mean(fit.loo.v1.m6$eval[,2,2]),
  mean(fit.loo.v1.m7$eval[,2,2]),
  mean(fit.loo.v1.m8$eval[,2,2]))

recall.v2 = c(mean(fit.loo.v2.m1$eval[,2,2]),
  mean(fit.loo.v2.m2$eval[,2,2]),
  mean(fit.loo.v2.m3$eval[,2,2]),
  mean(fit.loo.v2.m4$eval[,2,2]),
  mean(fit.loo.v2.m5$eval[,2,2]),
  mean(fit.loo.v2.m6$eval[,2,2]),
  mean(fit.loo.v2.m7$eval[,2,2]),
  mean(fit.loo.v2.m8$eval[,2,2]))


model = c("<(vd-2),<(vd-2)","1:nT,<(vd-2)","<pw,<(vd-2)",">=pw,<(vd-2)","<(vd-2),<pw","1:nT,<pw","<pw,<pw",">=pw,<pw")

recall.out=data.frame(model,recall.v1,recall.v2)
recall.out

###
### How do the epsilons change with different versions/models?
###
attributes(fit.loo.v1.m1)
epsilon.v1 = cbind(fit.loo.v1.m1$decide.star[,2],
              fit.loo.v1.m2$decide.star[,2],
              fit.loo.v1.m3$decide.star[,2],
              fit.loo.v1.m4$decide.star[,2],
              fit.loo.v1.m5$decide.star[,2],
              fit.loo.v1.m6$decide.star[,2],
              fit.loo.v1.m7$decide.star[,2],
              fit.loo.v1.m8$decide.star[,2])

epsilon.v2 = cbind(fit.loo.v2.m1$decide[,2],
              fit.loo.v2.m2$decide[,2],
              fit.loo.v2.m3$decide[,2],
              fit.loo.v2.m4$decide[,2],
              fit.loo.v2.m5$decide[,2],
              fit.loo.v2.m6$decide[,2],
              fit.loo.v2.m7$decide[,2],
              fit.loo.v2.m8$decide[,2])

epsilon.star.v1 = cbind(fit.loo.v1.m1$epsilon.star[,2],
               fit.loo.v1.m2$epsilon.star[,2],
               fit.loo.v1.m3$epsilon.star[,2],
               fit.loo.v1.m4$epsilon.star[,2],
               fit.loo.v1.m5$epsilon.star[,2],
               fit.loo.v1.m6$epsilon.star[,2],
               fit.loo.v1.m7$epsilon.star[,2],
               fit.loo.v1.m8$epsilon.star[,2])


epsilon.v1.df=data.frame(round(epsilon.v1,2))
epsilon.v2.df=data.frame(round(epsilon.v2,2))
epsilon.star.v1.df=data.frame(epsilon.star.v1)

model = c("<(vd-2),<(vd-2)","1:nT,<(vd-2)","<pw,<(vd-2)",">=pw,<(vd-2)","<(vd-2),<pw","1:nT,<pw","<pw,<pw",">=pw,<pw")

names(epsilon.v1.df)=model
names(epsilon.v2.df)=model
names(epsilon.star.v1.df)=model

epsilon.v1.df
epsilon.v2.df
epsilon.star.v1.df

recall.out=data.frame(model,recall.v1,recall.v2)
recall.out


##############################################################################################
###
### Final model fitting
### 
###
##############################################################################################



fit.recall.v1.m1 = training(d.dos,eps=fit.loo.v1.m1$decide.star[,2],pw=part.window,vd=vitdrop)# fit model v1 <(vd-2), test <vd-2
eval.recall.v1.m1=evaluate(alarm=fit.recall.v1.m1$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m1 = training(d.dos,eps=fit.loo.v2.m1$decide[,2],pw=part.window,vd=vitdrop)# fit model v1 <(vd-2), test <vd-2
eval.recall.v2.m1=evaluate(alarm=fit.recall.v2.m1$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

fit.recall.v1.m2 = training(d.dos,eps=fit.loo.v1.m2$decide.star[,2],pw=part.window,vd=vitdrop)# fit model #=1:nT, test <vd-2
eval.recall.v1.m2=evaluate(alarm=fit.recall.v1.m2$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m2 = training(d.dos,eps=fit.loo.v2.m2$decide[,2],pw=part.window,vd=vitdrop)# fit model #=1:nT, test <vd-2
eval.recall.v2.m2=evaluate(alarm=fit.recall.v2.m2$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

fit.recall.v1.m3 = training(d.dos,eps=fit.loo.v1.m3$decide.star[,2],pw=part.window,vd=vitdrop)# fit model <pw, test <vd-2
eval.recall.v1.m3=evaluate(alarm=fit.recall.v1.m3$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m3 = training(d.dos,eps=fit.loo.v2.m3$decide[,2],pw=part.window,vd=vitdrop)# fit model <pw, test <vd-2
eval.recall.v2.m3=evaluate(alarm=fit.recall.v2.m3$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

fit.recall.v1.m4 = training(d.dos,eps=fit.loo.v1.m4$decide.star[,2],pw=part.window,vd=vitdrop)# fit model >=pw:vd, test <vd-2
eval.recall.v1.m4=evaluate(alarm=fit.recall.v1.m4$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m4 = training(d.dos,eps=fit.loo.v2.m4$decide[,2],pw=part.window,vd=vitdrop)# fit model >=pw:vd, test <vd-2
eval.recall.v2.m4=evaluate(alarm=fit.recall.v2.m4$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)


fit.recall.v1.m5 = training(d.dos,eps=fit.loo.v1.m5$decide.star[,2],pw=part.window,vd=vitdrop)# fit model v1 <(vd-2), test <vd-2
eval.recall.v1.m5=evaluate(alarm=fit.recall.v1.m5$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m5 = training(d.dos,eps=fit.loo.v2.m5$decide[,2],pw=part.window,vd=vitdrop)# fit model v1 <(vd-2), test <vd-2
eval.recall.v2.m5=evaluate(alarm=fit.recall.v2.m5$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

fit.recall.v1.m6 = training(d.dos,eps=fit.loo.v1.m6$decide.star[,2],pw=part.window,vd=vitdrop)# fit model #=1:nT, test <vd-2
eval.recall.v1.m6=evaluate(alarm=fit.recall.v1.m6$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m6 = training(d.dos,eps=fit.loo.v2.m6$decide[,2],pw=part.window,vd=vitdrop)# fit model #=1:nT, test <vd-2
eval.recall.v2.m6=evaluate(alarm=fit.recall.v2.m6$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

fit.recall.v1.m7 = training(d.dos,eps=fit.loo.v1.m7$decide.star[,2],pw=part.window,vd=vitdrop)# fit model <pw, test <vd-2
eval.recall.v1.m7=evaluate(alarm=fit.recall.v1.m7$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m7 = training(d.dos,eps=fit.loo.v2.m7$decide[,2],pw=part.window,vd=vitdrop)# fit model <pw, test <vd-2
eval.recall.v2.m7=evaluate(alarm=fit.recall.v2.m7$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

fit.recall.v1.m8 = training(d.dos,eps=fit.loo.v1.m8$decide.star[,2],pw=part.window,vd=vitdrop)# fit model >=pw:vd, test <vd-2
eval.recall.v1.m8=evaluate(alarm=fit.recall.v1.m8$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2.m8 = training(d.dos,eps=fit.loo.v2.m8$decide[,2],pw=part.window,vd=vitdrop)# fit model >=pw:vd, test <vd-2
eval.recall.v2.m8=evaluate(alarm=fit.recall.v2.m8$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)



fit.recall = list(fit.recall.v1.m1,
                  fit.recall.v1.m2,
                  fit.recall.v1.m3,
                  fit.recall.v1.m4,
                  fit.recall.v1.m5,
                  fit.recall.v1.m6,
                  fit.recall.v1.m7,
                  fit.recall.v1.m8,
                  fit.recall.v2.m1,
                  fit.recall.v2.m2,
                  fit.recall.v2.m3,
                  fit.recall.v2.m4,
                  fit.recall.v2.m5,
                  fit.recall.v2.m6,
                  fit.recall.v2.m7,
                  fit.recall.v2.m8)


hit.recall=rep(list(),16)
output=function(x){
  hit.recall=array(NA,c(nInd.train,9))
  for(j in 1:nInd.train){
    hit.recall[j,1]=ifelse(!length(x$alarm[[j]]$julian[x$alarm[[j]]$julian <(vitdrop[j]-1)]),0.001,length(unique(x$alarm[[j]]$julian[x$alarm[[j]]$julian <(vitdrop[j]-1)])))
    hit.recall[j,2]=(vitdrop[j]-1)-126+1
    hit.recall[j,3]=hit.recall[j,1]/hit.recall[j,2]
    hit.recall[j,4]=sum(x$alarm[[j]]$julian <=vitdrop[j]+1 & x$alarm[[j]]$julian >=vitdrop[j]-1)>0
    hit.recall[j,5]=sum(x$alarm[[j]]$julian <=vitdrop[j]+2 & x$alarm[[j]]$julian >=vitdrop[j]-2)>0
    hit.recall[j,6]=sum(x$alarm[[j]]$julian <=vitdrop[j]+1 & x$alarm[[j]]$julian >=vitdrop[j]-1)
    hit.recall[j,7]=sum(x$alarm[[j]]$julian <=vitdrop[j]+2 & x$alarm[[j]]$julian >=vitdrop[j]-2)
    hit.recall[j,8]=sum(x$alarm[[j]]$julian ==vitdrop[j])
    hit.recall[j,9]=sum(x$alarm[[j]]$julian <(vitdrop[j]-1))
  }
  hit.recall
}

hit.recall = lapply(fit.recall,output)

day1.results = cbind(hit.recall[[1]][,4],
                     hit.recall[[2]][,4],
                     hit.recall[[3]][,4],
                     hit.recall[[4]][,4],
                     hit.recall[[5]][,4],
                     hit.recall[[6]][,4],
                     hit.recall[[7]][,4],
                     hit.recall[[8]][,4],
                     hit.recall[[9]][,4],
                     hit.recall[[10]][,4],
                     hit.recall[[11]][,4],
                     hit.recall[[12]][,4],
                     hit.recall[[13]][,4],
                     hit.recall[[14]][,4],
                     hit.recall[[15]][,4],
                     hit.recall[[16]][,4])

day2.results = cbind(hit.recall[[1]][,5],
                     hit.recall[[2]][,5],
                     hit.recall[[3]][,5],
                     hit.recall[[4]][,5],
                     hit.recall[[5]][,5],
                     hit.recall[[6]][,5],
                     hit.recall[[7]][,5],
                     hit.recall[[8]][,5],
                     hit.recall[[9]][,5],
                     hit.recall[[10]][,5],
                     hit.recall[[11]][,5],
                     hit.recall[[12]][,5],
                     hit.recall[[13]][,5],
                     hit.recall[[14]][,5],
                     hit.recall[[15]][,5],
                     hit.recall[[16]][,5])

day2.results
models=c("<(vd-2),<(vd-2)","1:nT,<(vd-2)","<pw,<(vd-2)",">=pw,<(vd-2)","<(vd-2),<pw","1:nT,<pw","<pw,<pw",">=pw,<pw","<(vd-2),<(vd-2)","1:nT,<(vd-2)","<pw,<(vd-2)",">=pw,<(vd-2)","<(vd-2),<pw","1:nT,<pw","<pw,<pw",">=pw,<pw")
results.variants = data.frame(cbind(c(rep("v1",8),rep("v2",8)),models,
                                    rbind(apply(hit.recall[[1]][,4:5],2,function(x){sum(x)/12}),
                                           apply(hit.recall[[2]][,4:5],2,function(x){sum(x)/12}),
                                           apply(hit.recall[[3]][,4:5],2,function(x){sum(x)/12}),
                                           apply(hit.recall[[4]][,4:5],2,function(x){sum(x)/12}),
                                           apply(hit.recall[[5]][,4:5],2,function(x){sum(x)/12}),
                                           apply(hit.recall[[6]][,4:5],2,function(x){sum(x)/12}),
                                           apply(hit.recall[[7]][,4:5],2,function(x){sum(x)/12}),
                                           apply(hit.recall[[8]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[9]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[10]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[11]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[12]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[14]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[14]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[15]][,4:5],2,function(x){sum(x)/12}),
                                          apply(hit.recall[[16]][,4:5],2,function(x){sum(x)/12}))),stringsAsFactors = FALSE)

names(results.variants) = c("Version","Model","Prop 1 Day","Prop 2 Day")
results.variants

results.fp = data.frame(cbind(hit.recall[[1]][,1],
                              hit.recall[[2]][,1],
                              hit.recall[[3]][,1],
                              hit.recall[[4]][,1],
                              hit.recall[[5]][,1],
                              hit.recall[[6]][,1],
                              hit.recall[[7]][,1],
                              hit.recall[[8]][,1],
                              hit.recall[[9]][,1],
                              hit.recall[[10]][,1],
                              hit.recall[[11]][,1],
                              hit.recall[[12]][,1],
                              hit.recall[[13]][,1],
                              hit.recall[[14]][,1],
                              hit.recall[[15]][,1],
                              hit.recall[[16]][,1]))
results.fp[results.fp==.001] = 0
results.fp
rep()
fp.summary=apply(results.fp,2,sum)
versions=c(rep("v1",8),rep("v2",8))
FP.df=data.frame(versions,model,fp.summary)

rownames(FP.df)=NULL

FP.df


###
### Alternative tables for evaluation (separated by training criteria)
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
hit.prec=matrix(NA,nr=nInd.train,nc = 9)
for(j in 1:nInd.train){
    hit.prec[j,1]=length(unique(fit.prec$alarm[[j]]$julian[fit.prec$alarm[[j]]$julian <(vitdrop[j]-1)]))
    hit.prec[j,2]=(vitdrop[j]-1)-126+1
    hit.prec[j,3]=hit.prec[j,1]/hit.prec[j,2]
    hit.prec[j,4]=sum(fit.prec$alarm[[j]]$julian <=vitdrop[j]+1 & fit.prec$alarm[[j]]$julian >=vitdrop[j]-1)>0
    hit.prec[j,5]=sum(fit.prec$alarm[[j]]$julian <=vitdrop[j]+2 & fit.prec$alarm[[j]]$julian >=vitdrop[j]-2)>0
    hit.prec[j,6]=sum(fit.prec$alarm[[j]]$julian <=vitdrop[j]+1 & fit.prec$alarm[[j]]$julian >=vitdrop[j]-1)
    hit.prec[j,7]=sum(fit.prec$alarm[[j]]$julian <=vitdrop[j]+2 & fit.prec$alarm[[j]]$julian >=vitdrop[j]-2)
    hit.prec[j,8]=sum(fit.prec$alarm[[j]]$julian ==vitdrop[j])
    hit.prec[j,9]=sum(fit.prec$alarm[[j]]$julian <(vitdrop[j]-1))
    }
hit.prec

write.csv(hit.prec,file="summary_hit_prec_v1_.csv")

hit.recall=matrix(NA,nr=n.vit,nc = 9)
for(j in 1:n.vit){
    hit.recall[j,1]=length(unique(fit.recall$alarm[[j]]$julian[fit.recall$alarm[[j]]$julian <(vitdrop[j]-1)]))
    hit.recall[j,2]=(vitdrop[j]-1)-126+1
    hit.recall[j,3]=hit.recall[j,1]/hit.recall[j,2]
    hit.recall[j,4]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+1 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-1)>0
    hit.recall[j,5]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+2 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-2)>0
    hit.recall[j,6]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+1 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-1)
    hit.recall[j,7]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+2 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-2)
    hit.recall[j,8]=sum(fit.recall$alarm[[j]]$julian ==vitdrop[j])
    hit.recall[j,9]=sum(fit.recall$alarm[[j]]$julian <(vitdrop[j]-1))
}
hit.recall

write.csv(hit.recall,file="summary_hit_recall_v2_m1.csv")

hit.recall=matrix(NA,nr=n.vit,nc = 9)
for(j in 1:n.vit){
  hit.recall[j,1]=length(unique(fit.recall$alarm[[j]]$julian[fit.recall$alarm[[j]]$julian <(vitdrop[j]-1)]))
  hit.recall[j,2]=(vitdrop[j]-1)-126+1
  hit.recall[j,3]=hit.recall[j,1]/hit.recall[j,2]
  hit.recall[j,4]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+1 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-1)>0
  hit.recall[j,5]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+2 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-2)>0
  hit.recall[j,6]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+1 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-1)
  hit.recall[j,7]=sum(fit.recall$alarm[[j]]$julian <=vitdrop[j]+2 & fit.recall$alarm[[j]]$julian >=vitdrop[j]-2)
  hit.recall[j,8]=sum(fit.recall$alarm[[j]]$julian ==vitdrop[j])
  hit.recall[j,9]=sum(fit.recall$alarm[[j]]$julian <(vitdrop[j]-1))
}
hit.recall

write.csv(hit.recall,file="summary_hit_recall_v2_m1.csv")




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

hit.recall

hit.F1=matrix(NA,nr=n.vit,nc = 9)
for(j in 1:n.vit){
    hit.F1[j,1]=length(unique(fit.F1$alarm[[j]]$julian[fit.F1$alarm[[j]]$julian <(vitdrop[j]-1)]))
    hit.F1[j,2]=(vitdrop[j]-1)-126+1
    hit.F1[j,3]=hit.F1[j,1]/hit.F1[j,2]
    hit.F1[j,4]=sum(fit.F1$alarm[[j]]$julian <=vitdrop[j]+1 & fit.F1$alarm[[j]]$julian >=vitdrop[j]-1)>0
    hit.F1[j,5]=sum(fit.F1$alarm[[j]]$julian <=vitdrop[j]+2 & fit.F1$alarm[[j]]$julian >=vitdrop[j]-2)>0
    hit.F1[j,6]=sum(fit.F1$alarm[[j]]$julian <=vitdrop[j]+1 & fit.F1$alarm[[j]]$julian >=vitdrop[j]-1)
    hit.F1[j,7]=sum(fit.F1$alarm[[j]]$julian <=vitdrop[j]+2 & fit.F1$alarm[[j]]$julian >=vitdrop[j]-2)
    hit.F1[j,8]=sum(fit.F1$alarm[[j]]$julian ==vitdrop[j])
    hit.F1[j,9]=sum(fit.F1$alarm[[j]]$julian <(vitdrop[j]-1))
}
hit.F1

write.csv(hit.F1,file="summary_hit_F1_v3.csv")


for(j in 1:10){
  beep(sound=5)
}



##############################################################################################
#SAVE working directory

save.image("anomaly_v3.Rdata")

#save
epsilon=epsilon.v2.df[,2]
save(epsilon,file='epsilon.Rdata')


##############################################################################################
##################################################################################################################################
###
### fitting altered versions of the threshold_detector.R function 
### to determine definition of data for quantile
###
##################################################################################################################################



detect.quant.test.m1=matrix(NA,nr=36,nc=16)
detect.quant.test.m2=matrix(NA,nr=36,nc=16)
detect.quant.test.m3=matrix(NA,nr=36,nc=16)
detect.quant.test.m4=matrix(NA,nr=36,nc=16)
detect.quant.test.m5=matrix(NA,nr=36,nc=16)
detect.quant.test.m6=matrix(NA,nr=36,nc=16)
detect.quant.test.m7=matrix(NA,nr=36,nc=16)
detect.quant.test.m8=matrix(NA,nr=36,nc=16)

for(k in 1:36){
  detect.quant.test.m1[k,]=fit.loo.m1$loo.fittest[[k]]$detect.quant
  detect.quant.test.m2[k,]=fit.loo.m2$loo.fittest[[k]]$detect.quant
  detect.quant.test.m3[k,]=fit.loo.m3$loo.fittest[[k]]$detect.quant
  detect.quant.test.m4[k,]=fit.loo.m4$loo.fittest[[k]]$detect.quant
  detect.quant.test.m5[k,]=fit.loo.m5$loo.fittest[[k]]$detect.quant
  detect.quant.test.m6[k,]=fit.loo.m6$loo.fittest[[k]]$detect.quant
  detect.quant.test.m7[k,]=fit.loo.m7$loo.fittest[[k]]$detect.quant
  detect.quant.test.m8[k,]=fit.loo.m8$loo.fittest[[k]]$detect.quant
}
detect.quant.test.m1
detect.quant.test.m2
detect.quant.test.m3
detect.quant.test.m4
detect.quant.test.m5
detect.quant.test.m6
detect.quant.test.m7
detect.quant.test.m8

fit.loo.m1$decide
nCovs=16
decide.versions = array(NA,c(nCovs,3,8))
decide.versions[,,1]=fit.loo.m1$decide
decide.versions[,,2]=fit.loo.m2$decide
decide.versions[,,3]=fit.loo.m3$decide
decide.versions[,,4]=fit.loo.m4$decide
decide.versions[,,5]=fit.loo.m5$decide
decide.versions[,,6]=fit.loo.m6$decide
decide.versions[,,7]=fit.loo.m7$decide
decide.versions[,,8]=fit.loo.m8$decide

decide.versions[,2,]




fit.recall.m1 = training(d.dos,eps=fit.loo.m1$decide[,2],pw=part.window,vd=vitdrop)# fit model v1 <(vd-2), test <vd-2
eval.recall.m1=evaluate(alarm=fit.recall.m1$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.m2 = training(d.dos,eps=fit.loo.m2$decide[,2],pw=part.window,vd=vitdrop)# fit model #=1:nT, test <vd-2
eval.recall.m2=evaluate(alarm=fit.recall.m2$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.m3 = training(d.dos,eps=fit.loo.m3$decide[,2],pw=part.window,vd=vitdrop)# fit model <pw, test <vd-2
eval.recall.m3=evaluate(alarm=fit.recall.m3$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.m4 = training(d.dos,eps=fit.loo.m4$decide[,2],pw=part.window,vd=vitdrop)# fit model >=pw:vd, test <vd-2
eval.recall.m4=evaluate(alarm=fit.recall.m4$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.m5 = training(d.dos,eps=fit.loo.m5$decide[,2],pw=part.window,vd=vitdrop)# fit model v1 <(vd-2), test <pw
eval.recall.m5=evaluate(alarm=fit.recall.m5$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.m6 = training(d.dos,eps=fit.loo.m6$decide[,2],pw=part.window,vd=vitdrop)# fit model #=1:nT, test <pw
eval.recall.m6=evaluate(alarm=fit.recall.m6$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.m7 = training(d.dos,eps=fit.loo.m7$decide[,2],pw=part.window,vd=vitdrop)# fit model <pw, test <pw
eval.recall.m7=evaluate(alarm=fit.recall.m7$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.m8 = training(d.dos,eps=fit.loo.m8$decide[,2],pw=part.window,vd=vitdrop)# fit model >=pw:vd, test <pw
eval.recall.m8=evaluate(alarm=fit.recall.m8$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)




fit.recall = list(fit.recall.m1,
                  fit.recall.m2,
                  fit.recall.m3,
                  fit.recall.m4,
                  fit.recall.m5,
                  fit.recall.m6,
                  fit.recall.m7,
                  fit.recall.m8)
hit.recall=rep(list(),8)
output=function(x){
  hit.recall=array(NA,c(nInd.train,9))
  for(j in 1:nInd.train){
    hit.recall[j,1]=ifelse(!length(x$alarm[[j]]$julian[x$alarm[[j]]$julian <(vitdrop[j]-1)]),0.001,length(unique(x$alarm[[j]]$julian[x$alarm[[j]]$julian <(vitdrop[j]-1)])))
    hit.recall[j,2]=(vitdrop[j]-1)-126+1
    hit.recall[j,3]=hit.recall[j,1]/hit.recall[j,2]
    hit.recall[j,4]=sum(x$alarm[[j]]$julian <=vitdrop[j]+1 & x$alarm[[j]]$julian >=vitdrop[j]-1)>0
    hit.recall[j,5]=sum(x$alarm[[j]]$julian <=vitdrop[j]+2 & x$alarm[[j]]$julian >=vitdrop[j]-2)>0
    hit.recall[j,6]=sum(x$alarm[[j]]$julian <=vitdrop[j]+1 & x$alarm[[j]]$julian >=vitdrop[j]-1)
    hit.recall[j,7]=sum(x$alarm[[j]]$julian <=vitdrop[j]+2 & x$alarm[[j]]$julian >=vitdrop[j]-2)
    hit.recall[j,8]=sum(x$alarm[[j]]$julian ==vitdrop[j])
    hit.recall[j,9]=sum(x$alarm[[j]]$julian <(vitdrop[j]-1))
  }
  hit.recall
}

hit.recall = lapply(fit.recall,output)

day1.results = cbind(hit.recall[[1]][,4],
                     hit.recall[[2]][,4],
                     hit.recall[[3]][,4],
                     hit.recall[[4]][,4],
                     hit.recall[[5]][,4],
                     hit.recall[[6]][,4],
                     hit.recall[[7]][,4],
                     hit.recall[[8]][,4])

day2.results = cbind(hit.recall[[1]][,5],
                     hit.recall[[2]][,5],
                     hit.recall[[3]][,5],
                     hit.recall[[4]][,5],
                     hit.recall[[5]][,5],
                     hit.recall[[6]][,5],
                     hit.recall[[7]][,5],
                     hit.recall[[8]][,5])

day2.results

results.variants = data.frame(cbind(c("<(vd-2),<(vd-2)","1:nT,<(vd-2)","<pw,<(vd-2)",">=pw,<(vd-2)","<(vd-2),<pw","1:nT,<pw","<pw,<pw",">=pw,<pw"),rbind(apply(hit.recall[[1]][,4:5],2,function(x){sum(x)/12}),
                                                                           apply(hit.recall[[2]][,4:5],2,function(x){sum(x)/12}),
                                                                           apply(hit.recall[[3]][,4:5],2,function(x){sum(x)/12}),
                                                                           apply(hit.recall[[4]][,4:5],2,function(x){sum(x)/12}),
                                                                           apply(hit.recall[[5]][,4:5],2,function(x){sum(x)/12}),
                                                                           apply(hit.recall[[6]][,4:5],2,function(x){sum(x)/12}),
                                                                           apply(hit.recall[[7]][,4:5],2,function(x){sum(x)/12}),
                                                                           apply(hit.recall[[8]][,4:5],2,function(x){sum(x)/12}))),stringsAsFactors = FALSE)

names(results.variants) = c("Model","Prop 1 Day","Prop 2 Day")
results.variants

results.fp = data.frame(cbind(hit.recall[[1]][,1],
                              hit.recall[[2]][,1],
                              hit.recall[[3]][,1],
                              hit.recall[[4]][,1],
                              hit.recall[[5]][,1],
                              hit.recall[[6]][,1],
                              hit.recall[[7]][,1],
                              hit.recall[[8]][,1]))
results.fp[results.fp==.001] = 0
names(results.fp)=c("<(vd-2),<(vd-2)","1:nT,<(vd-2)","<pw,<(vd-2)",">=pw,<(vd-2)","<(vd-2),<pw","1:nT,<pw","<pw,<pw",">=pw,<pw")
apply(results.fp,2,sum)


results.decide = data.frame(decide.versions[,2,])
names(results.decide)=c("<(vd-2),<(vd-2)","1:nT,<(vd-2)","<pw,<(vd-2)",">=pw,<(vd-2)","<(vd-2),<pw","1:nT,<pw","<pw,<pw",">=pw,<pw")
results.decide



detect.quant.test.m1=matrix(NA,nr=36,nc=16)
detect.quant.test.m2=matrix(NA,nr=36,nc=16)
detect.quant.test.m3=matrix(NA,nr=36,nc=16)
detect.quant.test.m4=matrix(NA,nr=36,nc=16)
detect.quant.test.m5=matrix(NA,nr=36,nc=16)
detect.quant.test.m6=matrix(NA,nr=36,nc=16)
detect.quant.test.m7=matrix(NA,nr=36,nc=16)
detect.quant.test.m8=matrix(NA,nr=36,nc=16)

for(k in 1:36){
  detect.quant.test.m1[k,]=fit.loo.m1$loo.fittest[[k]]$detect.quant
  detect.quant.test.m2[k,]=fit.loo.m2$loo.fittest[[k]]$detect.quant
  detect.quant.test.m3[k,]=fit.loo.m3$loo.fittest[[k]]$detect.quant
  detect.quant.test.m4[k,]=fit.loo.m4$loo.fittest[[k]]$detect.quant
  detect.quant.test.m5[k,]=fit.loo.m5$loo.fittest[[k]]$detect.quant
  detect.quant.test.m6[k,]=fit.loo.m6$loo.fittest[[k]]$detect.quant
  detect.quant.test.m7[k,]=fit.loo.m7$loo.fittest[[k]]$detect.quant
  detect.quant.test.m8[k,]=fit.loo.m8$loo.fittest[[k]]$detect.quant
}
detect.quant.test.m1
detect.quant.test.m2
detect.quant.test.m3
detect.quant.test.m4
detect.quant.test.m5
detect.quant.test.m6
detect.quant.test.m7
detect.quant.test.m8

length(fit.loo$loo.fittest[[1]]$detect.quant)

cbind(c(">=pw","<vd","nT","<=pw"),rbind(detect.quant.test[1,],detect.quant.test.m2[1,],detect.quant.test.m3[1,],detect.quant.test.m4[1,]))


epsilon=epsilon.v2.df[,2]

