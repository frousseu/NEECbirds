
### NEED TO GET OBJECTS FROM EBIRD!

library(RODBC)

db<-odbcConnectAccess("D:/ebird/Master ECSAS v 3.46.mdb")
obs<-sqlFetch(db,"tblSighting",as.is=TRUE) #?a plante et ne sait pas pourquoi
sp<-sqlFetch(db,"tblSpeciesInfo",as.is=TRUE) #?a plante et ne sait pas pourquoi
watch<-sqlFetch(db,"tblWatch",as.is=TRUE)
odbcClose(db)
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

x<-obs[,c("sp","date","month","year","lat","lon","nb")]

month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
names(month_comb)<-substr(names(month_comb),1,8)
x$season<-names(month_comb)[match(x$month,month_comb)]


### GET BIRD GROUPS
g<-fread("bird_groups.csv",sep=";",na.strings=c("","NA")) #ce fichier est sur mon github
g<-g[!is.na(group),]
g<-unique(g[,c("sp","group")])
x$group<-g$group[match(x$sp,g$sp)]



#################################
### grouping
#################################

group<-"seabirds_alcids"
sp<-unique(x$sp[x$group%in%group])
#sp<-"Dovekie"
month<-c("04","05","06","07")
m<-x$sp%in%sp & substr(x$date,6,7)%in%month
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
k<-kde(x=coordinates(xs),gridsize=c(500,500),compute.cont=TRUE,H=H1,w=xs$nb)
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
	season<-paste(month,collapse="")
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



















