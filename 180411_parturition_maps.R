###
### 1/31/2018 
### Alison Ketz
### Movement data preprocessing using adehabitatLT package
###


###
### Preliminaries
###

rm(list=ls())

library(e1071)
library(geosphere)
library(lubridate)
library(Hmisc)
library(zoo)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggmap)
library(adehabitatHR)
library(adehabitatLT)
library(maptools)
library(changepoint)
library(sp)
library(spatstat)#for "duplicate" function
library(readr)
library(RColorBrewer)


setwd("~/Documents/Parturition/180108_parturition")

# 
# par(mfrow=c(1,1))
# display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
#                    colorblindFriendly=TRUE)
# 

###
### Load VIT data
###

d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)[10]="lowtag"
d.vit$datedropped=as.character(d.vit$datedropped)
d.vit$juliandropped=yday(mdy_hms(d.vit$datedropped))
head(d.vit)
d.vit$datefawnfound=as.character(d.vit$datefawnfound)
d.vit$julianfawnfound=yday(mdy_hms(d.vit$datefawnfound))

n.vit=length(d.vit$lowtag)


#reorder vit data by lowtag/lowtag
d.vit=d.vit[order(d.vit$lowtag),]
 
#extract the individual ids
individs=d.vit$lowtag

#changing the 9th vit drop because it clearly dropped early
d.vit$juliandropped[9]=146 


###
### increased fixes parturition window
###

start=yday(mdy("03/22/2017")) # May 1 beginning of parturitionn period
end=yday(mdy("07/07/2017")) #end of parturition period

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


###
### double checking for duplicate locations
###

summary(duplicated(d))

#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_local))


###
### subset entire dataset to parturition window
###

d=d[d$julian>=start & d$julian <= end,]

#or keep all data, but exclude first couple hundred observations
#for plotting all data, not just window of parturition

# d=data.frame(d %>% group_by(lowtag) %>% slice(200:(n()-200)))
table(d$lowtag)

### Adding Vit data to main GPS dataframe

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
timediff <- diff(d$date_time_local)
d=d[-1,]
d$timediff <-round(as.numeric(abs(timediff)))
rm=which(d$timediff>10)
d=d[-rm,]

names(d)[1]="devname"
head(d)


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

tail(d)

#removing last day of parturition of 5004, outlier movements
rm5004 = (dim(d[d$lowtag==5004,])[1]-15):dim(d[d$lowtag==5004,])[1]
d=d[-rm5004,]

###
### Projections! 
###

# points from scratch
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
### The last movements of individ1 is crazy? running preditor?
### 
# plot(d.traj[1])
# 
# hist(d.traj[1], "dist", freq = TRUE)
# hist(d.traj[2], "dist", freq = TRUE)
# 
# plotltr(d.traj[1], "dt/3600/24")

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


#: the squared distance between the first relocation of the trajectory
# and the current relocation is often used to test some movements models
# (e.g. the correlated random walk, see the seminal paper of Kareiva and
# Shigesada, 1983).

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
#  Marsh and Jones 1988);

###
### mapping points using ggmap
###

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###
### Overall map of all individuals

# Make a bounding box for data
box <- make_bbox(longitude,latitude,data = d)
calc_zoom(box)
map=get_map(location=c(mean(d$longitude),mean(d$latitude)),source="google",zoom="auto",maptype="satellite",crop=box)
#Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=43.069755,-90.280971&zoom=14&size=640x640&scale=2&maptype=satellite&language=en-EN&sensor=false

pdf("Map_overall_Parturition.pdf")
ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=lowtag),data = d, size = 2)+
    scale_colour_brewer(palette = "Paired","ID")+xlab("Longitude")+ylab("Latitude")
dev.off()
###
### mapping points for each individual
###

n.map=length(individs)
d$dropped=as.factor(d$dropped)

for(i in 1:n.map){
    d.temp=d[d$lowtag==individs[i],]
    reflong=mean(d.temp$longitude)
    reflat=mean(d.temp$latitude)
    n.leg=dim(d.temp)[1]
    x=d.temp$longitude
    xend=c(d.temp$longitude[2:n.leg],d.temp$longitude[n.leg])
    y=d.temp$latitude
    yend=c(d.temp$latitude[2:n.leg],d.temp$latitude[n.leg])
    df.leg=data.frame(x,xend,y,yend)
    box <- make_bbox(longitude,latitude,data = d.temp)
    map=get_map(location=c(reflong,reflat),source="google",zoom=15,maptype="satellite",crop=box)
    pdf(file=paste("Map_",i,".pdf",sep=""))
    print(ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=dropped),data = d.temp, size = 2,alpha=.25)+
              scale_colour_manual(name = "", values = c(cbPalette[1],cbPalette[2]))+
              geom_point(aes(x=vitlongitude,y=vitlatitude),data=d.vit,size=4,alpha=1,color=cbPalette[3])+
              geom_point(aes(x = longitude, y = latitude),data = d.temp[d.temp$dropped=="1",], size = 2,alpha=1,color=cbPalette[2])+
              ggtitle(d.temp[1,1]))
    dev.off()
}


###
### subsetting the data to +- 3 days of the parturition event
###

d.subset=matrix(NA,nc=dim(d)[2]+1)
names(d.subset)=names(d)

for(j in 1:n.map){
    #subset data
    d.temp=d[d$lowtag==individs[j],]
    d.temp=d.temp[min(which(d.temp$julian==(d.vit$juliandropped[j]-1))):max(which(d.temp$julian==(d.vit$juliandropped[j]+1))),]
    d.temp$order=1:dim(d.temp)[1]
    # print with ggmap
    reflong=mean(d.temp$longitude)
    reflat=mean(d.temp$latitude)
    n.leg=dim(d.temp)[1]
    x=d.temp$longitude
    xend=c(d.temp$longitude[2:n.leg],d.temp$longitude[n.leg])
    y=d.temp$latitude
    yend=c(d.temp$latitude[2:n.leg],d.temp$latitude[n.leg])
    df.leg=data.frame(x,xend,y,yend)
    box <- make_bbox(longitude,latitude,data = d.temp)
    map=get_map(location=c(reflong,reflat),source="google",zoom=16,maptype="satellite",crop=box)
    pdf(file=paste("Map_dayof_",j,".pdf",sep=""))
    print(ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=dropped),data = d.temp, size = 2,alpha=.5)+
              scale_colour_manual(name = "", values = c(cbPalette[1],cbPalette[2]))+
              geom_path(aes(x = longitude, y = latitude),data = d.temp, size = 1,alpha=1,color=cbPalette[1])+
              geom_path(aes(x = longitude, y = latitude),data = d.temp[d.temp$dropped=="1",], size = 1,alpha=1,color=cbPalette[2])+
              geom_point(aes(x = longitude, y = latitude),data = d.temp[d.temp$dropped=="1",], size = 3,alpha=1,color=cbPalette[2])+
              geom_text(aes(x=longitude,y=latitude,label=order),data=d.temp[d.temp$dropped=="1",],size=2)+
              ggtitle(d.temp$lowtag[j]))
    
    dev.off()
    d.subset=rbind(d.subset,as.matrix(d.temp))
}

d.subset=d.subset[-1,]
names(d.subset)
names(d.subset)[]="order"
d.subset=data.frame(d.subset,stringsAsFactors = FALSE)
d.subset$date_time_local = as.POSIXct(d.subset$date_time_local)
d.subset$date_time_gmt = as.POSIXct(d.subset$date_time_gmt)
d.subset$coords.x1=as.numeric(d.subset$coords.x1)
d.subset$coords.x2=as.numeric(d.subset$coords.x2)


spdf.c=data.frame(d.subset[,20:21])
idsp=d.subset$lowtag
dtsp=d.subset$date_time_local
spdf=data.frame(idsp,dtsp)
coordinates(spdf)=spdf.c

#print with trajectory only$

d.sub.traj <- as.ltraj(coordinates(spdf),spdf$dtsp,id=spdf$idsp)
d.sub.traj
plot(d.sub.traj)
