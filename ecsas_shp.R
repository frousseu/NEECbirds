
### NEED TO GET OBJECTS FROM EBIRD!

library(RCurl)
library(RODBC)
library(plyr)
library(svMisc)
library(FRutils)
library(scales)

get_season<-function(x){
	s1<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
	s3<-unlist(list("12010203"=c("12","01","02","03"),"0405"=c("04","05"),"0607"=c("06","07"),"08091011"=c("08","09","10","11")))
	names(s1)<-substr(names(s1),1,nchar(names(s1))-1)
	names(s3)<-substr(names(s3),1,nchar(names(s3))-1)
	g<-grep("seabirds",x$group)	
	temp<-rep("not",nrow(x))
	if(any(g)){
		temp[g]<-"seabirds"
	}
	season<-ifelse(temp=="seabirds",names(s1)[match(x$month,s1)],names(s3)[match(x$month,s3)])
	season
}

ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
laea<-"+proj=laea +lat_0=50 +lon_0=-65"

### get ECSAS data
#db<-odbcConnectAccess("D:/ebird/Master ECSAS v 3.51.mdb")
#obs<-sqlFetch(db,"tblSighting",as.is=TRUE)
#sp<-sqlFetch(db,"tblSpeciesInfo",as.is=TRUE)
#watch<-sqlFetch(db,"tblWatch",as.is=TRUE)
#odbcClose(db)
load("D:/ebird/odbcECSAS.RData") #load session instead of running odbcConnect on 32bit all the time
rm(kde2pol)


### get SOMEC data
#db<-odbcConnectAccess2007("D:/ebird/SOMEC_CopieMOD31012017_clean.accdb")
load("D:/ebird/SOMECexport2011andover.RData") #load session instead of running on 32bit all the time
addQC$WatchID<-paste0("qc",addQC$WatchID)
addQC$Date<-substr(addQC$"Start Date",1,10)
watchQC<-addQC[!duplicated(addQC$WatchID),c("WatchID","Date","LongStart","LatStart")]
obsQC<-addQC[!addQC$Alpha%in%c("",NA),]


### join ECSAS and SOMEC data
watch<-join(watch,watchQC,type="full")
obs<-join(obs,obsQC,type="full")

obs$Alpha<-NA
obs$Alpha<-ifelse(is.na(obs$Alpha),sp$Alpha[match(obs$SpecInfoID,sp$SpecInfoID)],obs$Alpha)
obs$sp<-sp$English[match(obs$Alpha,sp$Alpha)]
obs$Class<-sp$Class[match(obs$Alpha,sp$Alpha)]
obs<-obs[obs$Class=="Bird",]
m<-match(obs$WatchID,watch$WatchID)
obs$date<-substr(watch$Date[m],1,10)
obs$lat<-watch$LatStart[m]
obs$lon<-watch$LongStart[m]
obs$lat<-ifelse(!is.na(obs$ObsLat),obs$ObsLat,obs$lat)
obs$lon<-ifelse(!is.na(obs$ObsLong),obs$ObsLong,obs$lon)
obs$lat<-ifelse(obs$lat<0,-obs$lat,obs$lat)
obs$lon<-ifelse(obs$lon>0,-obs$lon,obs$lon)

obs<-obs[obs$lon>=(-80) & obs$lon<=30,]
obs<-obs[obs$lat>=(40) & obs$lat<=75,]

obs$month<-substr(obs$date,6,7)
obs$year<-substr(obs$date,1,4)
obs<-obs[obs$year>="2000",]
obs$nb<-obs$Count
obs$nb<-ifelse(is.na(obs$nb),1,obs$nb)

#month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
#names(month_comb)<-substr(names(month_comb),1,8)
#obs$season<-names(month_comb)[match(obs$month,month_comb)]

### GET BIRD GROUPS
g<-getURL("https://raw.githubusercontent.com/frousseu/NEECbirds/master/bird_groups.csv") # Ce fichier est sur mon github
g<-read.csv(text=g,header=TRUE,stringsAsFactors=FALSE)
g<-g[!is.na(g$group_kernels),]
g<-unique(g[,c("sp","group_kernels")])
obs$group<-g$group_kernels[match(obs$sp,g$sp)]
obs$season<-get_season(obs)


### keep only on water obs for waterfowl
nofly<-c("waterfowl_diving","waterfowl_dabbling")
w<-which(obs$group%in%nofly & obs$FlySwim=="F")
obs<-obs[-w,]

x<-obs[,c("sp","group","season","date","month","year","lat","lon","nb")]


#################################################################    
### get seasonal effort and evaluate at each observation location
#################################################################



watch<-watch[substr(watch$Date,1,4)>="2000",]
watch$month<-substr(watch$Date,6,7)
#watch$season<-names(month_comb)[match(watch$month,month_comb)]
ys<-watch
xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon,x$lat)),ncol=2),proj4string=CRS(ll),data=data.frame(id=row.names(x),month=x$month,stringsAsFactors=FALSE))
xs<-spTransform(xs,CRS(prj))
ys$lat<-ys$LatStart
ys$lon<-ys$LongStart
ys<-ys[ys$lon>=(-100) & ys$lon<=30,]
ys<-ys[ys$lat>=(40) & ys$lat<=75,]
ys<-SpatialPointsDataFrame(matrix(as.numeric(c(ys$lon,ys$lat)),ncol=2),proj4string=CRS(ll),data=data.frame(month=ys$month))
ys<-spTransform(ys,CRS(prj))

weights<-unlist(lapply(unique(get_season(x)),function(i){
	mm<-unlist(strsplit(gsub("(.{2})", "\\1 ",i)," "))
	xxs<-xs[xs$month%in%mm,]	
	yys<-ys[ys$month%in%mm,]	
	H1<-matrix(c(200000000,0,0,200000000),nrow=2)  
	ke<-kde(x=coordinates(yys),binned=TRUE,eval.points=coordinates(xxs),compute.cont=FALSE,H=H1) # binned accélère BEAUCOUP le calcul et ne semb<le pas changer grand chose
	res<-ke$estimate
	names(res)<-xxs$id
	res
}))

weights<-weights[match(xs$id,names(weights))]
x$we<-1-(weights/max(weights))



#################################
### grouping
#################################

### transform weights
f<-function(x){log(x+1)}
#f<-function(x){x}

case<-ddply(x,.(group,season),nrow)
case<-case[!case$group%in%c("passerines","raptors",NA) & case$V1>=100,]

for(j in seq_len(nrow(case))){
	
	group<-case$group[j]
	season<-case$season[j]
	#sp<-"Dovekie"
	m<-x$group%in%group & x$season%in%season
	xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll),data=x[m,])
	xs<-spTransform(xs,CRS(prj))
	
	
	######################################
	### KERNELS KDE
	######################################
	
	#H<-Hpi.diag(coordinates(xs))
	#H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
	#H1<-H*matrix(c(0.02,0,0,0.02),nrow=2) 
	H1<-matrix(c(200000000,0,0,200000000),nrow=2) 
	### GET KERNELS POLYGON
	k<-kde(x=coordinates(xs),binned=FALSE,gridsize=c(500,500),compute.cont=TRUE,H=H1,w=f(xs$nb)*xs$we)
	kp<-list()
	perc<-c(25,50,75,95)
	percw<-c("very high","high","medium","low")
	trans<-c(0.8,0.6,0.4,0.2)
	cols_kern<-c("darkred","red","orange","yellow")
	kp<-kde2pol(k,levels=perc,proj4string=proj4string(xs),cut=FALSE) # extract polygons
 kp$group<-group
 kp$season<-season
 kp$site<-NA
 kp$sp<-NA
	
 png(paste0("D:/ebird/kernels/",paste(group,season,sep="_"),"_ecsas.png"),width=12,height=8,units="in",res=500,pointsize=14)
	
 ### PLOT POLYGONS
	par(mar=c(1,0,0,0))
	plot(xs,col="white")
	plot(kp,add=TRUE,col=alpha(cols_kern,0.6),border=NA)
	plot(land,col="grey95",border="grey75",add=TRUE)
	legend("bottomright",fill=alpha(cols_kern,0.6),legend=paste(perc,"%"),border=NA,cex=1,bg="grey90",box.lwd=NA,inset=c(0.05,0.1),title="Kernel Contours")
	info<-paste0("group: ",group,"\nseason: ",season,"\ndata source: ECSAS","\ncoast buffer:")
	mtext(info,side=3,line=-4,font=2,adj=0.05)
	
	dev.off()
	
	kp$data_type<-"ecsas"
	kp$risk_level<-percw
	
	### write shapefile
	writeOGR(kp,dsn="D:/ebird/kernels",layer=paste(group,season,"ecsas",sep="_"),driver="ESRI Shapefile",overwrite_layer=TRUE)
	
	progress(j,nrow(case))
	
}

















