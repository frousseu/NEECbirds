library(data.table)
library(sp)
library(rgdal)
library(rgeos)
library(scales)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(raster)
library(ks)
library(FRutils)
library(svMisc)
library(RODBC)
library(readxl)
library(plyr)

d<-list()
month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
names(month_comb)<-substr(names(month_comb),1,8)


######################################
### NOGA
######################################

x<-list.files("D:/Télémétrie/Donnees_GPS")
x<-x[-grep("ali",x)]
l<-sapply(x,function(i){
	path<-paste0("D:/Télémétrie/Donnees_GPS/",i) 
	sh<-excel_sheets(path)
})
l<-data.frame(f=rep(names(l),sapply(l,length)),s=unlist(l))
l<-apply(l,1,function(i){
	ans<-read_excel(path=paste0("D:/Télémétrie/Donnees_GPS/",i[1]),sheet=i[2])
	ans<-ans[,1:5]
	ans$id<-i[2]
	ans
})
x<-do.call("rbind",l)
x<-x[x$Comportement%in%c("pêche","dérive"),]
x<-as.data.frame(x)
x<-x[!is.na(x$Latitude),]
x$year<-substr(x$Date,1,4)
x$y<-as.numeric(x$Latitude)
x$x<-as.numeric(x$Longitude)
x$id2<-1
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
x<-x[,c("id","Date","Latitude","Longitude")]
names(x)<-c("id","date","lat","lon")
x$sp<-"NOGA"
n<-ddply(x@data,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$season<-names(month_comb)[match(x$month,month_comb)]
d[["NOGA"]]<-x


##############################################
### RTLO
##############################################

x<-read.csv("D:/Télémétrie/RTLO_EEZ.txt",stringsAsFactors=FALSE)
x$x<-x$longitud
x$y<-x$latitude
x$id<-x$animal
x$id2<-1
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
x<-x[,c("animal","date","latitude","longitud")]
names(x)<-c("id","date","lat","lon")
x$date<-substr(strptime(x$date,"%m/%d/%Y %H:%M:%S"),1,10)
x$sp<-"RTLO"
n<-ddply(x@data,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$season<-names(month_comb)[match(x$month,month_comb)]
d[["RTLO"]]<-x



###########################################
### RAZO
###########################################

#db<-odbcConnectAccess2007("D:/Télémétrie/ID_General_RAZO.accdb")
#d2015<-sqlFetch(db,"Donnees_2015",stringsAsFactors=FALSE)
#d2016<-sqlFetch(db,"Donnees_2016",stringsAsFactors=FALSE)
#odbcClose(db)
d2015<-as.data.frame(fread("D:/Télémétrie/Donnees_2015_RAZO.txt",stringsAsFactors=FALSE,na.strings=c("","NA")))
d2016<-as.data.frame(fread("D:/Télémétrie/Donnees_2016_RAZO.txt",stringsAsFactors=FALSE,na.strings=c("","NA")))
names(d2015)[which(names(d2015)=="Longtitude")]<-"Longitude"
x<-full_join(d2015,d2016)
x<-x[!is.na(x$Latitude),]
x$x<-as.numeric(x$Longitude)
x$y<-as.numeric(x$Latitude)
x$id<-x$LoggerID
x$id2<-1
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
x<-x[x$id%in%names(table(x$id)[table(x$id)>2]),]
x<-x[,c("LoggerID","date","Latitude","Longitude")]
names(x)<-c("id","date","lat","lon")
x$date<-substr(strptime(x$date,"%m/%d/%Y %H:%M:%S"),1,10)
x$sp<-"RAZO"
n<-ddply(x@data,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$season<-names(month_comb)[match(x$month,month_comb)]
d[["RAZO"]]<-x


#############################################
### GSGO
#############################################

x<-as.data.frame(fread("D:/Télémétrie/data2012_01_17_OIES.txt",stringsAsFactors=FALSE))
x<-x[x$typ=="gps" & x$xx%in%0:5 & x$good,]
x<-x[x$lat>45 & x$lat<50,]
x<-x[x$lon>(-77),]
x$x<-x$lon
x$y<-x$lat
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
x<-x[,c("band","date","lat","lon")]
names(x)<-c("id","date","lat","lon")
x$sp<-"GSGO"
n<-ddply(x@data,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$season<-names(month_comb)[match(x$month,month_comb)]
d[["GSGO"]]<-x


#####################################
### BUILD SEASONAL SUBSETS
#####################################

season<-unique(names(month_comb))
d<-unlist(lapply(d,function(i){
     res<-lapply(season,function(j){
    	  ans<-i[i$season==j,]
    	  if(nrow(ans)==0){
    	    "empty"
    	  }else{
    	    ans
    	  }
     })
     names(res)<-season
     res
}))
names(d)<-gsub("\\.","_",names(d))
d<-d[sapply(d,function(i){!identical(i,"empty")})]




####################################################
### GET KERNELS
####################################################


for(j in seq_along(d)){
	 x<-d[[j]]
	 group<-strsplit(names(d)[j],"_")[[1]][1]
	 season<-strsplit(names(d)[j],"_")[[1]][2]
	 #H<-Hpi.diag(coordinates(x))
	 #H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
	 #H1<-H*matrix(c(0.02,0,0,0.02),nrow=2) 
	 H1<-matrix(c(20000000,0,0,20000000),nrow=2) 
  k<-kde(x=coordinates(x),binned=TRUE,bgridsize=c(500,500),compute.cont=TRUE,H=H1,w=x$we)
  perc<-c(25,50,75,95)
  percw<-c("very high","high","medium","low")
  trans<-c(0.9,0.7,0.5,0.3)
  cols_kern<-c("darkred","red","orange","yellow")
  kp<-kde2pol(k,levels=perc,proj4string=prj) # extract polygons
  kp$group<-names(d)[j]
  kp$season<-season
  
  
  ####################################################################
  ### PLOT PNG
  ####################################################################
  
  png(paste0("D:/ebird/kernels/",paste(group,season,sep="_"),"_tracking.png"),width=12,height=8,units="in",res=600,pointsize=14)
  
  #par(mar=c(8.5,6,0,0),mfrow=c(1,1))
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(bbox2pol(x),border="white")
  plot(land,col="grey95",border="grey75",add=TRUE,lwd=0.5)
  
  ### KERNS
  plot(kp,add=TRUE,col=alpha(cols_kern,0.6),border=NA)
  #plot(gIntersection(coast,kp,byid=TRUE),add=TRUE,col=alpha(cols_kern,0.6),border=NA)
  info<-paste0("group: ",group,"\nseason: ",season,"\ndata source: tracking")
  mtext(info,side=3,line=-4,font=2,adj=0.05)
  legend("topleft",title="Kernel Contours (risk)",fill=alpha(cols_kern,0.6),legend=paste(perc,"%","(",percw,")"),border=NA,cex=1,,bg="grey90",box.lwd=NA,inset=c(0.05,0.2))
  
  dev.off()
  
  ##########################################################
  ### WRITE SHAPEFILE
  ###########################################################
  
  kp$data_type<-"tracking"
  kp$risk_level<-percw
  
  writeOGR(kp,dsn="D:/ebird/kernels",layer=paste(group,season,"tracking",sep="_"),driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  progress(j,length(d))
  
  
}



















