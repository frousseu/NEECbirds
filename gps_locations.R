library(readxl)
library(adehabitatHR)
library(ks)
library(sp)
library(rgdal)
library(rgeos)
library(splancs)
library(svMisc)
library(leaflet)
library(magrittr)
library(RODBC)
library(FRutils)
library(dplyr)
library(hexbin)
library(scales)

#save.image("W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Code_R/gps_locations.RData")
#load("W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Code_R/gps_locations.RData") #contient maintenant le preleaflet2

load("C:/Users/User/Documents/SCF2016_FR/UrgencesAviR/Dev/urgencesapp/www/urgencesapp.RData")


colo<-function(x=0:1,rescale=TRUE,col=c("white","yellow","orange","red","darkred")){
  if(rescale){
    x<-x/max(x)
  }
  palette<-colorRamp(col)(x)
  palette<-rgb(palette[,1],palette[,2],palette[,3],maxColorValue=256)
  palette
}

kde2pol<-function(k,perc="5%",proj=NULL){
  co<-with(k,contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate,levels=cont[perc]))
  poly<-lapply(co,function(i){
    x<-data.frame(i$x,i$y)
    x<-rbind(x,x[1,]) 
    x<-SpatialPolygons(list(Polygons(list(Polygon(x)),ID=1)))
  })
  poly<-lapply(seq_along(poly),function(i){
    spChFIDs(poly[[i]], as.character(i))
  })
  poly<-do.call("rbind",poly)
  proj4string(poly)<-CRS(proj)
  poly
}


projected<-"+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"


na<-readOGR(dsn="C:/Users/User/Documents/SCF2016_FR/shapefiles",layer="ne_10m_admin_1_states_provinces",encoding="UTF-8")
na<-na[na$geonunit%in%c("Greenland","Canada"),]
na<-gIntersection(na,bbox2pol(na,ex=-1),byid=TRUE)
na<-spTransform(na,CRS(projected))



### FOUS
x<-list.files("C:/Users/User/Documents/SCF2016_FR/Télémétrie/Donnees_GPS")
x<-x[-grep("ali",x)]
l<-sapply(x,function(i){
  path<-paste0("C:/Users/User/Documents/SCF2016_FR/Télémétrie/Donnees_GPS/",i) 
  sh<-excel_sheets(path)
})
l<-data.frame(f=rep(names(l),sapply(l,length)),s=unlist(l))
l<-apply(l,1,function(i){
  ans<-read_excel(path=paste0("C:/Users/User/Documents/SCF2016_FR/Télémétrie/Donnees_GPS/",i[1]),sheet=i[2])
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
noga<-x



### RTLO
x<-read.csv("C:/Users/User/Documents/SCF2016_FR/Télémétrie/RTLO_EEZ.txt",stringsAsFactors=FALSE)
x$x<-x$longitud
x$y<-x$latitude
x$id<-x$animal
x$id2<-1
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
rtlo<-x
#range(substr(strptime(rtlo$date,"%m/%d/%Y %H:%M:%S")),6,10)

### RAZO
db<-odbcConnectAccess2007("C:/Users/User/Documents/SCF2016_FR/Télémétrie/ID_General_RAZO.accdb")
d2015<-sqlFetch(db,"Donnees_2015",stringsAsFactors=FALSE)
d2016<-sqlFetch(db,"Donnees_2016",stringsAsFactors=FALSE)
odbcClose(db)
names(d2015)[which(names(d2015)=="Longtitude")]<-"Longitude"
x<-full_join(d2015,d2016)
x<-x[!is.na(x$Latitude),]
x$x<-x$Longitude
x$y<-x$Latitude
x$id<-x$LoggerID
x$id2<-1
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
x<-x[x$id%in%names(table(x$id)[table(x$id)>2]),]
razo<-x

### GSGO
x<-read.table("C:/Users/User/Documents/OIES/Data/data2012_01_17.txt",header=TRUE,stringsAsFactors=FALSE)
x<-x[x$typ=="gps" & x$xx%in%0:5 & x$good,]
x<-x[x$lat>45 & x$lat<50,]
x<-x[x$lon>(-77),]
coordinates(x)<-~lon+lat
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
gsgo<-x







### CHOOSE SPECIES

x<-razo

### BUILD GRID

n<-150 # number of divisions wished in the convex hull bbox
region<-gConvexHull(x)
b<-bbox(region)
s<-round(abs(b[1,2]-b[1,1])/n,0)
g<-GridTopology(c(b[1,1],b[2,1]),c(s,s),c(ceiling((b[1,2]-b[1,1])/s),ceiling((b[2,2]-b[2,1])/s)))
g<-SpatialGrid(g)
grid<-g
grid<-as(grid,"SpatialPolygons")
grid2<-spsample(grid,type="hexagonal",cellsize=s)
grid<-HexPoints2SpatialPolygons(grid2)
grid<-SpatialPolygonsDataFrame(grid,data=data.frame(id=1:length(grid)),match.ID=FALSE)
proj4string(grid)<-CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80")
o<-over(grid,region)
grid<-grid[!is.na(o),]
g1<-!apply(gContains(gUnaryUnion(na),grid,byid=TRUE),1,any)
grid<-grid[g1,]


### BUILD KERNESL

kmaster<-kde(x=coordinates(x),compute.cont=TRUE)
#r<-raster(kmaster)
#plot(r)
ind<-names(table(x$id)[table(x$id)>10])  ### on subset ce qui a plus de 10 localisations
hr<-vector(mode="list",length=length(ind))
names(hr)<-ind
for(i in seq_along(ind)){
 #i<-40
 H<-Hpi.diag(coordinates(x[x$id==ind[i],]))
 H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic)
 H<-H*matrix(c(0.5,0,0,0.5),nrow=2) 
 #H<-Hpi.diag(coordinates(x[x$id==ind[i],]))
 k<-kde(x=coordinates(x[x$id==ind[i],]),compute.cont=TRUE,H=H)
 kp<-kde2pol(k,perc="5%",proj=proj4string(x))
 hr[[i]]<-kp
 #plot(region,border="white")
 #plot(municip.shp.proj,add=TRUE,col=ifelse(is.na(municip.shp.proj$MUS_NM_MUN),"lightblue","grey"))
 #plot(k, display="filled.contour2", cont=seq(5,95,by=5),asp=1,add=TRUE)
 #plot(k,cont=seq(5,5,by=5),asp=1,add=TRUE)
 #plot(kp,border="blue",add=TRUE)
 #plot(x[x$id==ind[i],],pch=16,cex=0.5,add=TRUE,col=alpha("black",0.99))
 progress(i,length(ind))
}

lo<-NULL
for(i in seq_along(hr)){
  y<-over(grid,hr[[i]])
  lo<-c(lo,grid$id[!is.na(y)])
  progress(i,length(hr))
}
tab<-table(lo)
grid$p<-tab[match(grid$id,names(tab))]/length(hr)
grid$p<-ifelse(is.na(grid$p),0,grid$p)
#par(mar=c(0,0,0,6))
#plot(x,col="white")
#plot(municip.shp.proj,add=TRUE,col=ifelse(is.na(municip.shp.proj$MUS_NM_MUN),"lightblue2","grey55"),border=NA)
#plot(grid[grid$p>0,],col=colo(grid$p[grid$p>0]),border=NA,add=TRUE)
#legend("topright",inset=c(-0.0,0.25),legend=paste(round(100*seq(0,max(grid$p),length.out=6),0),"%"),fill=colo(seq(0,max(grid$p),length.out=6)),border=NA,cex=1.5,bg="white",xpd=TRUE)
#plot(municip.shp.proj,add=TRUE)



par(mar=c(0,0,0,0))
plot(region,border="white")
for(i in seq_along(hr)){
	plot(hr[[i]],col=alpha("red",0.1),border=NA,add=TRUE)
	progress(i,length(hr))
}


### relation between nb locations and home range surface
nbloc<-sapply(ind,function(i){sum(i==x$id)})
surf<-sapply(ind,function(i){gArea(hr[[i]])/(1000*1000)})
plot(nbloc,surf)


rtlo_shp<-grid[grid$p>0,]
rtlo_shp$p<-as.numeric(rtlo_shp$p)
rtlo_shp_ll<-spTransform(rtlo_shp,CRS("+proj=longlat +datum=WGS84"))

#writeOGR(rtlo_shp_ll,dsn="W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Code_R/shiny/www",layer="RTLO_shp",driver="ESRI Shapefile",overwrite_layer=TRUE)

r<-range(rtlo_shp_ll$p)
se<-seq(r[1],r[2],length.out=6)
leaflet() %>%
  addProviderTiles("Esri.WorldImagery",options = providerTileOptions(noWrap = TRUE)) %>%
  addPolygons(data=rtlo_shp_ll,stroke=FALSE,fillColor=colo(rtlo_shp_ll$p),popup=rtlo_shp_ll$p,weight=0,fillOpacity=0.4,opacity=0) %>%
  addLegend(values=se,color=colo(se),labels=format(100*se,digits=1,nsmall=1),opacity=0.4,position="bottomright",title="% of individual kernel homeranges touching cell")







### get grid using kernels from adehabitat
#res<-kernelUD(x[,"id2"],h="href",same4all=TRUE,grid=200,extent=1)
#grid<-res[[1]]@grid
#grid<-SpatialGrid(grid,proj4string=CRS(proj4string(x)))
#grid<-as(grid,"SpatialPixels")
#grid<-as(grid,"SpatialPolygons")
#grid<-SpatialPolygonsDataFrame(grid,data=data.frame(id=names(grid),stringsAsFactors=F),match.ID=F)
#hr<-getverticeshr(res,95)























