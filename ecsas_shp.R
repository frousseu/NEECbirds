
### NEED TO GET OBJECTS FROM EBIRD!

library(RODBC)
library(plyr)
library(svMisc)

db<-odbcConnectAccess("D:/ebird/Master ECSAS v 3.46.mdb")
obs<-sqlFetch(db,"tblSighting",as.is=TRUE)
sp<-sqlFetch(db,"tblSpeciesInfo",as.is=TRUE)
watch<-sqlFetch(db,"tblWatch",as.is=TRUE)
odbcClose(db)

load("D:/ebird/odbcECSAS.RData") #load session instead of running on 32bit all the time

obs$Alpha<-sp$Alpha[match(obs$SpecInfoID,sp$SpecInfoID)]
obs$sp<-sp$English[match(obs$SpecInfoID,sp$SpecInfoID)]
obs$Class<-sp$Class[match(obs$SpecInfoID,sp$SpecInfoID)]
obs<-obs[obs$Class=="Bird",]
m<-match(obs$WatchID,watch$WatchID)
obs$date<-substr(watch$Date[m],1,10)
obs$lat<-watch$LatStart[m]
obs$lon<-watch$LongStart[m]
obs$lat<-ifelse(!is.na(obs$ObsLat),obs$ObsLat,obs$lat)
obs$lon<-ifelse(!is.na(obs$ObsLong),obs$ObsLong,obs$lon)
obs$lat<-ifelse(obs$lat<0,-obs$lat,obs$lat)
obs$lon<-ifelse(obs$lon>0,-obs$lon,obs$lon)

obs<-obs[obs$lon>=(-100) & obs$lon<=30,]
obs<-obs[obs$lat>=(40) & obs$lat<=75,]

obs$month<-substr(obs$date,6,7)
obs$year<-substr(obs$date,1,4)
obs<-obs[obs$year>="2000",]
obs$nb<-obs$Count

month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
names(month_comb)<-substr(names(month_comb),1,8)
obs$season<-names(month_comb)[match(obs$month,month_comb)]

### GET BIRD GROUPS
g<-fread("bird_groups.csv",sep=";",na.strings=c("","NA")) #ce fichier est sur mon github
g<-g[!is.na(group),]
g<-unique(g[,c("sp","group")])
obs$group<-g$group[match(obs$sp,g$sp)]

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
watch$season<-names(month_comb)[match(watch$month,month_comb)]
ys<-watch
xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon,x$lat)),ncol=2),proj4string=CRS(ll),data=data.frame(id=row.names(x),season=x$season,stringsAsFactors=FALSE))
xs<-spTransform(xs,CRS(prj))
ys$lat<-ys$LatStart
ys$lon<-ys$LongStart
ys<-ys[ys$lon>=(-100) & ys$lon<=30,]
ys<-ys[ys$lat>=(40) & ys$lat<=75,]
ys<-SpatialPointsDataFrame(matrix(as.numeric(c(ys$lon,ys$lat)),ncol=2),proj4string=CRS(ll),data=data.frame(season=ys$season))
ys<-spTransform(ys,CRS(prj))

weights<-unlist(lapply(unique(names(month_comb)),function(i){
	xxs<-xs[xs$season==i,]	
	yys<-ys[ys$season==i,]	
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
	sp<-"Dovekie"
	m<-x$sp%in%sp & x$season%in%season
	xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll),data=x[m,])
	xs<-spTransform(xs,CRS(prj))
	
	
	######################################
	### KERNELS KDE
	######################################
	
	#H<-Hpi.diag(coordinates(xs))
	#H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
	#H1<-H*matrix(c(0.02,0,0,0.02),nrow=2) 
	H1<-matrix(c(200000000,0,0,200000000),nrow=2) 
	
	### GET KERNELS POLYGONS
	k<-kde(x=coordinates(xs),binned=FALSE,gridsize=c(500,500),compute.cont=TRUE,H=H1,w=f(xs$nb)*xs$we)
	kp<-list()
	perc<-c(25,50,75,95)
	percw<-c("very high","high","medium","low")
	trans<-c(0.8,0.6,0.4,0.2)
	cols_kern<-c("darkred","red","orange","yellow")
	for(i in seq_along(perc)){
		kp[[i]]<-kde2pol(k,perc=paste0(100-perc[i],"%"),proj=proj4string(xs)) # extract polygons
	}
	for(i in rev(seq_along(kp))){
		if(i==1){ #make sure a single polygon for each contour
			kp[[i]]<-gUnaryUnion(kp[[i]])
		}else{	
			kp[[i]]<-gSymdifference(kp[[i]],kp[[i-1]],byid=FALSE) # keep non overlapping parts
		}
		id<-paste0("k",perc[i])
		res<-SpatialPolygonsDataFrame(kp[[i]],data=data.frame(id=id,group=group,season=season,stringsAsFactors=FALSE),match.ID=FALSE)
		kp[[i]]<-spChFIDs(res,id) # give unique ID
	}
	kp<-do.call("rbind",kp)
	
	### PLOT POLYGONS
	par(mar=c(1,0,0,0))
	plot(xs,col="white")
	plot(kp,add=TRUE,col=alpha(cols_kern,0.6),border=NA)
	plot(land,col="grey95",border="grey75",add=TRUE)
	legend("bottomright",fill=alpha(cols_kern,0.6),legend=paste(perc,"%"),border=NA,cex=1,bg="grey90",box.lwd=NA,inset=c(0.05,0.1),title="Kernel Contours")
	
	### write shapefile
	writeOGR(kp,dsn="D:/ebird/kernels",layer=paste(group,season,"ecsas",sep="_"),driver="ESRI Shapefile",overwrite_layer=TRUE)
	
	progress(j,nrow(case))
	
}

















