library(ff)
library(data.table)
library(sp)
library(rgdal)
library(rgeos)
library(quantreg)
library(scales)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(spatstat)
library(raster)
library(OpenStreetMap)
library(ks)
library(FRutils)
library(yarrr)

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

load("D:/ebird/birdgroupings.rda")
g<-birdgroupings
rm(birdgroupings)
g$wetland<-apply(g[,c("Habitat1","Habitat2","Habitat3","Habitat4")],1,function(i){any(i=="we")})
g<-g[,c("English_Name","Species_gr2")]
names(g)<-c("sp","group")

#con<-file("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt")
#open(con)
#d<-read.table(con,skip=5000000,nrow=1,sep="\t")

co<-readOGR(dsn="D:/ebird/ebd_CA_relAug-2016.txt",layer="ne_10m_admin_1_states_provinces",encoding="UTF-8")

na<-readOGR(dsn="D:/ebird/ebd_CA_relAug-2016.txt",layer="ne_10m_admin_1_states_provinces",encoding="UTF-8")
na<-na[na$admin%in%c("Canada","United States of America"),]
na<-na[na$name%in%c("Quebec","Nova Scotia","Vermont","New York","New Hampshire","QuÃ©bec","Prince Edward Island","New Brunswick","Newfoundland and Labrador","Maine"),]

land<-spTransform(na,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))


g1<-gBuffer(gUnaryUnion(land),width=-5000)
g2<-gBuffer(gUnaryUnion(land),width=5000)
coast<-gDifference(g2,g1)
coastw<-gDifference(coast,land)
plot(land)
plot(coastw,add=TRUE,col="blue",border=NA)




ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=22 +datum=NAD83 +ellps=GRS80"
laea<-"+proj=laea +lat_0=50 +lon_0=-65"

keep<-c("COMMON NAME","TAXONOMIC ORDER","STATE","LATITUDE","LONGITUDE","OBSERVATION DATE","OBSERVATION COUNT","CATEGORY","SAMPLING EVENT IDENTIFIER")
keepn<-c("sp","taxo","state","lat","lon","date","nb","category","sample")

# 29 409 415 rows in there
e<-fread("D:/ebird/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",encoding = "UTF-8", na.strings = "",header = TRUE, sep = "\t",colClasses = rep("character", 44),verbose = TRUE,select=keep,nrows=29409415)
names(e)<-keepn[match(names(e),keep)]


e<-merge(e,g,all=TRUE)
e<-e[e$category%in%c("species"),] #on perd des formes et qques taxons/formes
e$month<-substr(e$date,6,7)
#e<-as.data.table(e)


month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
names(month_comb)<-substr(names(month_comb),1,8)

e$season<-names(month_comb)[match(e$month,month_comb)]

res<-e[,.(nbobs=.N),by=.(group,state,season)]
res<-res[!is.na(group)]
res$state_season<-paste0(res$state,res$season)
res<-res[res$group%in%c("Auks and allies","Cormorants","Dabbling","Diving","Geese","Gulls","Jaegers","Loons and Grebes","Petrels and allies","Shorebirds","Terns"),]
#res<-res[grep("New|Edw|Scot",res$state),]
res<-spread(res[,c("group","state_season","nbobs")],state_season,nbobs,fill=0)

options(scipen=20)
par(mar=c(18,7,6,2))
barplot(apply(res[,-1],2,identity),border="white",las=2,col=brewer.pal(nrow(res), "Paired"))
legend("topright",legend=rev(res$group),fill=rev(brewer.pal(nrow(res), "Paired")))
mtext("Nb of records",2,line=5)
mtext("Nb of ebird records per group and season",3,line=3)








#fwrite(e[e$STATE=="Quebec",],"D:/ebird/ebd_CA_relAug-2016.txt/ebirdQC.csv",sep="\t")

e<-fread("D:/ebird/ebd_CA_relAug-2016.txt/ebirdQC.csv",encoding = "UTF-8", na.strings = "",header = TRUE, sep = "\t",verbose = TRUE)

names(e)<-keepn[match(names(e),keep)]

x<-e[state%in%c("Quebec","Nova Scotia","New Brunswick","Prince Edward Island","Newfoundland and Labrador"),]
x<-x[category%in%c("species"),]

x$taxo<-as.numeric(x$taxo)
temp<-unique(x[,c("sp","taxo")])
temp<-temp[order(temp$taxo),]
x$nb<-as.numeric(x$nb)
x$nb<-ifelse(is.na(x$nb),1,x$nb)
x$lon<-as.numeric(x$lon)
x$lat<-as.numeric(x$lat)
x<-x[lat<=52,] #on garde ce qui est en bas


#################################
### grouping
#################################

sp<-"Semipalmated Sandpiper"
group<-"Auks and allies"
month<-c("07","08","09","10")
m<-x$sp%in%sp & substr(x$date,6,7)%in%month
xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll),data=x[m,])
xs<-spTransform(xs,CRS(proj4string(land)))
o<-over(xs,coast)
xs<-xs[!is.na(o),]
grid<-FRutils:::hexgrid(xs,width=NULL,n=100,convex=TRUE)
#boxplot(nb~o$id)

nb<-xs$nb
lat<-coordinates(xs)[,2]
lon<-coordinates(xs)[,1]

xx<-as.data.table(xs@data)
o<-over(xs,grid)
xx$id<-o$id
xx[,n:=.N,by=id] #pass by reference no assignment
#ids<-unique(xx$id[which(xx$nb>50 | xx$n>100)])
#pirateplot(nb~id,data=xx[xx$id%in%ids,],quant=0.9)

### keep coastal part of grid
o<-over(grid,coast)
grid<-grid[!is.na(o),]

### add zeros
o<-apply(over(grid,xs),1,function(j){all(is.na(j))})
addlat<-jitter(rep(coordinates(grid[o,])[,2],times=10),amount=5000)
addlon<-jitter(rep(coordinates(grid[o,])[,1],times=10),amount=5000)
addnb<-rep(0,length(addlat))
lat<-c(lat,addlat)
lon<-c(lon,addlon)
nb<-c(nb,addnb)

#windows()


#####################################
### RQSS
#####################################

fit<-rqss(nb~qss(cbind(lon,lat),lambda=0.005),tau=0.80)

g<-gBuffer(gConvexHull(xs),width=-1000)
o1<-over(SpatialPoints(coordinates(grid),CRS(proj4string(grid))),g)
o2<-apply(over(grid,xs),1,function(j){all(is.na(j))})
w<-!is.na(o1) & !is.na(o2) 
X<-coordinates(grid)[w,1]
Y<-coordinates(grid)[w,2]
p<-predict(fit,data.frame(lon=X,lat=Y,stringsAsFactors=FALSE))[,1]
p<-ifelse(p<0,0,p)

cols<-rev(c("white","yellow","red","darkred"))
#trans<-function(x,min=0.05){ans<-sqrt(x)/max(sqrt(x));ans<-ifelse(ans<min,min,ans);ans}

layout(matrix(c(1,2),ncol=1))
par(mar=c(1,0,0,0))
plot(xs,col="white")
plot(land,col="grey95",border="grey75",add=TRUE)
plot(grid[w,],col=colo.scale(sqrt(p),cols),border=NA,add=TRUE)
leg<-(seq(sqrt(1),sqrt(max(p)),length.out=12))^2
col<-tail(colo.scale(sqrt(c(p,leg)),cols),length(leg))
leg<-paste0(c("\u2264",rep("",length(leg)-1)),round(leg,0))
legend("right",legend=rev(leg),fill=rev(col),cex=1.5,border=NA,bg="grey90",box.lwd=NA,inset=c(0.05,0))
plot(land,border="grey75",add=TRUE)

#plot(g1,border="green",add=TRUE)
#plot(g2,border="red",add=TRUE)

######################################
### KERNELS KDE
######################################

H<-Hpi.diag(coordinates(xs))
H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
#H1<-H*matrix(c(0.5,0,0,0.5),nrow=2) 
H1<-H*matrix(c(0.1,0,0,0.1),nrow=2) 
#H1<-matrix(c(50000000,0,0,50000000),nrow=2) 

### GET POLYGONS
k<-kde(x=coordinates(xs),compute.cont=TRUE,H=H1,w=nb)
kp<-list()
perc<-c(25,50,75,95)
trans<-c(0.8,0.6,0.4,0.2)
cols<-c("darkred","red","orange","yellow")
for(i in seq_along(perc)){
  kp[[i]]<-kde2pol(k,perc=paste0(100-perc[i],"%"),proj=proj4string(grid)) # extract polygons
}
for(i in rev(seq_along(kp)[-1])){
	kp[[i]]<-gSymdifference(kp[[i]],kp[[i-1]]) # keep non overlapping parts
}

### PLOT POLYGONS
par(mar=c(1,0,0,0))
plot(xs,col="white")
plot(land,col="grey95",border="grey75",add=TRUE)
for(i in rev(seq_along(perc))){
	plot(kp[[i]],add=TRUE,col=alpha("red",trans[i]),border=NA)
}
#plot(xs,add=TRUE,pch=1,col=alpha("black",0.25),cex=15*nb/max(nb))
legend("right",fill=alpha("red",trans),legend=paste(perc,"%"),border=NA,cex=1.5,,bg="grey90",box.lwd=NA,inset=c(0.05,0))
#
#

### BARPLOT ABUNDANCE
par(mar=c(2,2,0,0),new=TRUE)
h<-hist(xs$nb,breaks=seq(0,max(xs$nb),by=20),plot=FALSE)
h<-hist(xs$nb,breaks=seq(0,max(xs$nb),by=20),col=alpha("red",0.5),border=ifelse(h$counts,"red",NA),xaxt="n",yaxt="n",)
axis(1)
axis(2)
#
#
#


####################################################################
### FINAL KERNEL OUTPUT
####################################################################

### plot the look of the output
par(mar=c(4,4,0,0),mfrow=c(1,1))
plot(xs,col="white")
plot(land,col="grey95",border="grey75",add=TRUE)
b<-bbox2pol(proj4string=proj4string(coastw))
plot(gIntersection(coastw,b,byid=TRUE),col=alpha("blue",0.15),border=NA,add=TRUE)
for(i in rev(seq_along(perc))){
  plot(kp[[i]],add=TRUE,col=cols[i],border=NA)
	 #plot(gIntersection(coastw,kp[[i]],byid=TRUE),add=TRUE,col=cols[i],border=NA)
}
#plot(xs,add=TRUE,pch=1,col=alpha("black",0.25),cex=15*nb/max(nb))

### get an histogram of the size of each record
h<-lapply(kp,function(i){
  o<-over(xs,i)
  res<-xs$nb[!is.na(o)]
  brks<-c(0,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000)
  brks<-brks[brks<=max(xs$nb)]
  res<-cut(res,breaks=brks)
  table(res)
})

### add the histogram to the plot
par(new=TRUE)
barplot(do.call("rbind",h),col=alpha(cols,0.35),beside=TRUE,border=NA,xlab="Classe d'abondance",ylab="Nombre de mentions")

### LEAFLET

kpc<-lapply(kp,function(i){gIntersection(coastw,i,byid=TRUE)})


leaflet() %>%
	addProviderTiles("Esri.WorldImagery",options = providerTileOptions(noWrap = TRUE)) %>%
	#addPolygons(data=spTransform(kpc[[1]],ll),stroke=FALSE,fillColor=cols[1],weight=0,fillOpacity=0.7,opacity=0) %>% 
 #addPolygons(data=spTransform(kpc[[2]],ll),stroke=FALSE,fillColor=cols[2],weight=0,fillOpacity=0.7,opacity=0) %>% 
 #addPolygons(data=spTransform(kpc[[3]],ll),stroke=FALSE,fillColor=cols[3],weight=0,fillOpacity=0.7,opacity=0) %>% 
 #addPolygons(data=spTransform(kpc[[4]],ll),stroke=FALSE,fillColor=cols[4],weight=0,fillOpacity=0.7,opacity=0) %>% 
	#addPolygons(data=na,stroke=FALSE,fillColor="white",weight=0,fillOpacity=0.5,opacity=0) %>% 
addPolylines(data=co)
	#addCircleMarkers(data=spTransform(xs[1:1000,],ll),stroke=FALSE,fillColor=cols[4],weight=0,fillOpacity=0.7,opacity=0)




#####################################
### KERNELS spatstat
#####################################

b<-bbox(xs)
hppp<-ppp(lon,lat,xrange=c(b[1,1],b[1,2]),yrange=c(b[2,1],b[2,2]),marks=nb)
spatstat.options(npixel=c(1000,1000))
h<-density(hppp,sigma=0.005,weights=marks(hppp),positive=FALSE)
m<-t(h$v)
val<-quantile(m,0.997)
m[]<-ifelse(m<val,NA,val)
dat1=list()
dat1$x=seq(h$xcol[1],by=h$xstep,len=h$dim[1])
dat1$y=seq(h$yrow[1],by=h$ystep,len=h$dim[2])
dat1$z=m
r<-raster(dat1)
#plot(map.osm) 
#plot(r,bty="n",axes=FALSE,box=FALSE,xpd=FALSE,col="red")
#plot(fleuve,add=TRUE)#,border=gray(0.5,0.3))
#text(par("usr")[1]+4,par("usr")[4]-0,paste(c("nb mentions:",nrow(d2),"range:",range(d2$Abundance)),collapse=" "),xpd=TRUE,adj=c(0,1))
pr<-rasterToPolygons(r,dissolve=TRUE)
pr<-spChFIDs(pr,sp)
res<-spTransform(pr,CRS(proj4string(grid)))
plot(pr,add=TRUE,col=alpha("darkgreen",0.75),border=NA)
#plot(spTransform(pr,osm()),col=alpha("red",0.5),border=NA)
















































