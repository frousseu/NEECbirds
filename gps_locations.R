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
mrc.shp <- readOGR(dsn = "W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Data_output/Quebec",layer = "mrc_s")
municip.shp <- readOGR(dsn = "W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Data_output/Quebec",layer= "munic_s") 
na <- readOGR(dsn = "W:/PIU/Projet_BanqueDeDonnees/SHAPEFILES/PoliticalBoundaries_Shapefiles/NA_PoliticalDivisions/data/bound_p",layer="boundary_p_v2")
usa<-na[sort(unique(c((1:nrow(na))[na$NAME%in%c("Quebec / Québec","Nova Scotia / Nouvelle-????cosse","Nunavut","Newfoundland and Labrador / Terre-Neuve-et-Labrador","New Brunswick / Nouveau-Brunswick","Prince Edward Island / ?Zle-du-Prince-????douard","Vermont","New Hampshire","New York","Maine","Massachusetts","Rhode Island","Delaware")],grep("Scotia|Prince",na$NAME)))) ,]
#usa<-na[sort(unique(c((1:nrow(na))[na$NAME%in%c("Nova Scotia / Nouvelle-????cosse","Nunavut","Newfoundland and Labrador / Terre-Neuve-et-Labrador","New Brunswick / Nouveau-Brunswick","Prince Edward Island / ?Zle-du-Prince-????douard","Vermont","New Hampshire","New York","Maine","Massachusetts","Rhode Island","Delaware")],grep("Scotia|Prince",na$NAME)))) ,]

usa<-spTransform(usa,CRS(projected))
municip.shp.proj<-spTransform(municip.shp,CRS(projected))
mrc.shp.proj<-spTransform(mrc.shp,CRS(projected))
east<-gUnaryUnion(municip.shp.proj[is.na(municip.shp.proj$MUS_NM_MUN),],id=rep(1,nrow(municip.shp.proj[is.na(municip.shp.proj$MUS_NM_MUN),])))

na<-readOGR(dsn="C:/Users/User/Documents/SCF2016_FR/shapefiles",layer="ne_10m_admin_1_states_provinces",encoding="UTF-8")
na<-na[na$geonunit%in%c("Greenland","Canada"),]
na<-gIntersection(na,bbox2pol(na,ex=-1),byid=TRUE)
na<-spTransform(na,CRS(projected))


### build region of interest

windows()
par(mar=c(0,0,0,0))
plot(na)
region <- getpoly()
region<-region[c(1:nrow(region),1),]
region<-SpatialPolygons(list(Polygons(list(Polygon(region)),ID=1)),proj4string=CRS(proj4string(municip.shp.proj)))
b<-bbox(region)


s<-10000
#g<-GridTopology(c(b[1,1],b[2,1]),c(s,s),c(ceiling((b[1,2]-b[1,1])/s),ceiling((b[2,2]-b[2,1])/s)))
g<-GridTopology(c(b[1,1],b[2,1]),c(s,s),c(ceiling((b[1,2]-b[1,1])/s),ceiling((b[2,2]-b[2,1])/s)))
g<-SpatialGrid(g)
grid<-g
grid<-as(grid,"SpatialPolygons")
grid<-SpatialPolygonsDataFrame(grid,data=data.frame(id=1:length(grid)),match.ID=FALSE)
proj4string(grid)<-CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80")
o<-over(grid,region)
grid<-grid[!is.na(o),]


g1<-!apply(gIntersects(gUnaryUnion(na),grid,byid=TRUE),1,any)

plot(grid)
plot(grid[g1,],col="blue",add=TRUE)


#################
#g1<-!apply(gIntersects(gUnaryUnion(usa),grid,byid=TRUE),1,any)



#################

#plot(grid,border="red")
#o2<-apply(over(grid,usa),1,function(i){any(!is.na(i))})
#grid<-grid[o2,]

#plot(grid,border="blue")
#o3<-over(grid,municip.shp.proj[is.na(municip.shp.proj$MUS_NM_MUN),])
#plot(grid,border="green")

#g<-gCrosses(usa,grid,byid=TRUE)
#g<-gContains(gUnaryUnion(usa),grid,byid=TRUE)
#grid<-grid[!apply(g,1,any),]



#grid<-grid[!is.na(o),]

#g<-gIntersection(grid,usa,unaryUnion_if_byid_false=FALSE)


#plot(grid,border="white")
#plot(usa,add=TRUE)
#plot(municip.shp.proj,add=TRUE,border="blue")
#plot(grid,border="red",add=TRUE)

#o<-apply(over(grid,usa),1,function(i){any(!is.na(i))})
#grid2<-grid[o,]
#plot(grid2,add=TRUE,col="green")
#g2<-gIntersection(grid2,usa)









#fleuve<-municip.shp.proj[is.na(municip.shp.proj$MUS_NM_MUN),]
#o<-!apply(over(grid,fleuve),1,function(i){all(is.na(i))})
#g<-gIntersection(grid[o,],fleuve,byid=TRUE)


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
#o<-over(x,region)
#x<-x[!is.na(o),]


### RTLO
x<-read.csv("C:/Users/User/Documents/SCF2016_FR/Télémétrie/RTLO_EEZ.txt",stringsAsFactors=FALSE)
x$x<-x$longitud
x$y<-x$latitude
x$id<-x$animal
x$id2<-1
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
o<-over(x,region)
x<-x[!is.na(o),]
rtlo<-x
#range(substr(strptime(rtlo$date,"%m/%d/%Y %H:%M:%S")),6,10)

### RAZO
db<-odbcConnectAccess2007("C:/Users/User/Documents/SCF2016_FR/Télémétrie/ID_General_RAZO.accdb")
d2015<-sqlFetch(db,"Donnees_2015",stringsAsFactors=FALSE)
d2016<-sqlFetch(db,"Donnees_2016",stringsAsFactors=FALSE)
odbcClose(db)
names(d2015)[which(names(d2015)=="Longtitude")]<-"Longitude"
x<-full_join(d2015,d2016)
#x<-read.csv("W:/PIU/DONNÉES DE RÉPARTITION DES OISEAUX/GPS2015.csv",stringsAsFactors=FALSE)
x<-x[!is.na(x$Latitude),]
x$x<-x$Longitude
x$y<-x$Latitude
x$id<-x$LoggerID
x$id2<-1
coordinates(x)<-~x+y
proj4string(x)<-CRS("+proj=longlat +datum=NAD83")
x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
o<-over(x,region)
x<-x[!is.na(o),]
x<-x[x$id%in%names(table(x$id)[table(x$id)>2]),]
razo<-x








### CHOOSE SPECIES

x<-noga 



### BUILD GRID

n<-100
region<-gConvexHull(x)
b<-bbox(region)
s<-round(abs(b[1,2]-b[1,1])/n,0)
g<-GridTopology(c(b[1,1],b[2,1]),c(s,s),c(ceiling((b[1,2]-b[1,1])/s),ceiling((b[2,2]-b[2,1])/s)))
g<-SpatialGrid(g)
grid<-g
grid<-as(grid,"SpatialPolygons")
grid<-SpatialPolygonsDataFrame(grid,data=data.frame(id=1:length(grid)),match.ID=FALSE)
proj4string(grid)<-CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80")
o<-over(grid,region)
grid<-grid[!is.na(o),]
g1<-!apply(gContains(gUnaryUnion(na),grid,byid=TRUE),1,any)
grid<-grid[g1,]
plot(grid)






### BUILD KERNESL

kmaster<-kde(x=coordinates(x),compute.cont=TRUE)
#r<-raster(kmaster)
#plot(r)
ind<-unique(x$id)
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
par(mar=c(0,0,0,6))
plot(x,col="white")
plot(municip.shp.proj,add=TRUE,col=ifelse(is.na(municip.shp.proj$MUS_NM_MUN),"lightblue2","grey55"),border=NA)
plot(grid[grid$p>0,],col=colo(grid$p[grid$p>0]),border=NA,add=TRUE)
legend("topright",inset=c(-0.0,0.25),legend=paste(round(100*seq(0,max(grid$p),length.out=6),0),"%"),fill=colo(seq(0,max(grid$p),length.out=6)),border=NA,cex=1.5,bg="white",xpd=TRUE)
plot(municip.shp.proj,add=TRUE)




rtlo_shp<-grid[grid$p>0,]
rtlo_shp$p<-as.numeric(rtlo_shp$p)
rtlo_shp_ll<-spTransform(rtlo_shp,CRS("+proj=longlat +datum=WGS84"))

#writeOGR(rtlo_shp_ll,dsn="W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Code_R/shiny/www",layer="RTLO_shp",driver="ESRI Shapefile",overwrite_layer=TRUE)

leaflet() %>%
  addProviderTiles("Stamen.TonerLite",options = providerTileOptions(noWrap = TRUE)) %>%
	addTiles(group="Base") %>%
  addPolygons(data=rtlo_shp_ll,stroke=FALSE,fillColor=colo(rtlo_shp_ll$p),popup=rtlo_shp_ll$p,weight=0,fillOpacity=0.6,opacity=0)
  #colorNumeric("Reds",domain=c(0,1))







### get grid using kernels from adehabitat
#res<-kernelUD(x[,"id2"],h="href",same4all=TRUE,grid=200,extent=1)
#grid<-res[[1]]@grid
#grid<-SpatialGrid(grid,proj4string=CRS(proj4string(x)))
#grid<-as(grid,"SpatialPixels")
#grid<-as(grid,"SpatialPolygons")
#grid<-SpatialPolygonsDataFrame(grid,data=data.frame(id=names(grid),stringsAsFactors=F),match.ID=F)
#hr<-getverticeshr(res,95)























