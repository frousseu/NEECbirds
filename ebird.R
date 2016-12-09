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

ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
laea<-"+proj=laea +lat_0=50 +lon_0=-65"

### GET BIRD GROUPS
g<-fread("bird_groups.csv",sep=";",na.strings=c("","NA")) #ce fichier est sur mon github
g<-g[!is.na(group),]
g<-unique(g[,c("sp","group")])


#co<-readOGR(dsn="D:/ebird/ebd_CA_relAug-2016.txt",layer="ne_10m_admin_1_states_provinces",encoding="UTF-8")
#co<-readOGR(dsn="D:/ebird",layer="coastlines_z1",encoding="UTF-8")

na<-readOGR(dsn="D:/ebird/ebd_CA_relAug-2016.txt",layer="ne_10m_admin_1_states_provinces",encoding="UTF-8")
land<-na[na$admin%in%c("Canada","United States of America","Greenland"),]
land<-land[land$name%in%c("Manitoba","Virginia","Delaware","Ontario","Nunavut","Quebec","Nova Scotia","Vermont","New York","New Hampshire","QuÃ©bec","Prince Edward Island","New Brunswick","Newfoundland and Labrador","Maine","Massachusetts","Maryland","Pennsylvania","Connecticut","New Jersey","Rhode Island","Ohio","Kentuchy","West Virginia","Greenland") | land$admin%in%"Greenland",]
#na<-na[na$name%in%c("British Columbia"),]

qcmaritimes<-na[na$admin%in%c("Canada","United States of America"),]
qcmaritimes<-qcmaritimes[qcmaritimes$name%in%c("Ontario","Quebec","Nova Scotia","Vermont","New York","New Hampshire","QuÃ©bec","Prince Edward Island","New Brunswick","Newfoundland and Labrador","Maine"),]

land<-spTransform(land,CRS(prj))
qcmaritimes<-spTransform(qcmaritimes,CRS(prj)) # smaller subset to make the coast cause it's faster

g1<-gBuffer(gUnaryUnion(gBuffer(qcmaritimes,width=1)),width=-5000)
g2<-gBuffer(gUnaryUnion(gBuffer(qcmaritimes,width=1)),width=5000)
coast<-gDifference(g2,g1)
coastw<-gDifference(coast,qcmaritimes)


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

#e<-fread("D:/ebird/ebd_CA_relAug-2016.txt/ebirdQC.csv",encoding = "UTF-8", na.strings = "",header = TRUE, sep = "\t",verbose = TRUE)
#names(e)<-keepn[match(names(e),keep)]

x<-e[state%in%c("Quebec","Nova Scotia","New Brunswick","Prince Edward Island","Newfoundland and Labrador"),]
#x<-e[state%in%c("British Columbia"),]


x$taxo<-as.numeric(x$taxo)
temp<-unique(x[,c("sp","taxo")])
temp<-temp[order(temp$taxo),]
x$nb<-as.numeric(x$nb)
x$nb<-ifelse(is.na(x$nb),1,x$nb)
x$lon<-as.numeric(x$lon)
x$lat<-as.numeric(x$lat)
x<-x[lat<=52,] #on garde ce qui est en bas
x<-x[lon>=(-74),] #on garde ce qui est à partir de mtl


#################################
### EFFORT
#################################

xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon,x$lat)),ncol=2),proj4string=CRS(ll),data=data.frame(id=row.names(x),season=x$season,stringsAsFactors=FALSE))
xs<-spTransform(xs,CRS(prj))

y<-unique(x[,c("sample","lat","lon")])
ys<-SpatialPoints(matrix(as.numeric(c(y$lon,y$lat)),ncol=2),proj4string=CRS(ll))
ys<-spTransform(ys,CRS(prj))
o<-over(ys,coast)
ys<-ys[!is.na(o),]

H1<-matrix(c(6000000,0,0,6000000),nrow=2) 

### GET POLYGONS
#k<-kde(x=coordinates(ys),binned=TRUE,bgridsize=c(200,200),compute.cont=TRUE,H=H1)

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



kp<-list()
perc<-c(seq(5,95,by=10),99)
trans<-rev(seq(0.15,1,length.out=length(perc)))
col_eff<-c(colo.scale(perc[-length(perc)],c("blue","violet","magenta","yellow")),alpha("blue",0.5))
for(i in seq_along(perc)){
	kp[[i]]<-kde2pol(k,perc=paste0(100-perc[i],"%"),proj=proj4string(xs)) # extract polygons
}
for(i in rev(seq_along(kp)[-1])){
	kp[[i]]<-gSymdifference(kp[[i]],kp[[i-1]],byid=FALSE) # keep non overlapping parts
}

### PLOT POLYGONS
par(mar=c(1,0,0,0))
plot(bbox2pol(kp[[length(kp)]]),col="white",border="white")
plot(land,border="grey75",add=TRUE)
for(i in seq_along(kp)){
	#plot(kp[[i]],add=TRUE,col=alpha("red",trans[i]),border=NA)
	plot(gIntersection(coast,kp[[i]]),add=TRUE,col=col_eff[i],border=NA)
}
legend("bottomright",fill=col_eff,legend=paste(perc,"%"),border=NA,cex=1,box.lwd=NA,inset=c(0.05,0.1),title="Kernel Contours")



#################################
### grouping
#################################

group<-"seabirds_alcids"
#sp<-"Glaucous-winged Gull"
sp<-unique(x$sp[x$group%in%group])
month<-c("04","05","06","07")
season<-paste(month,collapse="")
m<-x$sp%in%sp & substr(x$date,6,7)%in%month
xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll),data=x[m,])
xs<-spTransform(xs,CRS(prj))
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
windows(w=21,h=12)

tau<-0.90
fit<-rqss(nb~qss(cbind(lon,lat),lambda=0.05),tau=tau)

g<-gBuffer(gConvexHull(xs),width=-1000)
o1<-over(SpatialPoints(coordinates(grid),CRS(proj4string(grid))),g)
o2<-apply(over(grid,xs),1,function(j){all(is.na(j))})
w<-!is.na(o1) & !is.na(o2) 
X<-coordinates(grid)[w,1]
Y<-coordinates(grid)[w,2]
p<-predict(fit,data.frame(lon=X,lat=Y,stringsAsFactors=FALSE))[,1]
p<-ifelse(p<0,0,p)

cols_gam<-rev(c("lightgreen","yellow","orange","red","darkred"))
#trans<-function(x,min=0.05){ans<-sqrt(x)/max(sqrt(x));ans<-ifelse(ans<min,min,ans);ans}

f<-function(x){log(x+5)} #transformation to better appreciate variation in group sizes
finv<-function(x){exp(x)-5}


layout(matrix(c(1,2),ncol=1))
par(mar=c(1,0,0,0))
plot(xs,col="white")
plot(land,col="grey95",border="grey75",add=TRUE)
plot(grid[w,],col=colo.scale(f(p),cols_gam),border=NA,add=TRUE)
leg<-finv((seq(f(1),f(max(p)),length.out=12)))
col<-tail(colo.scale(f(c(p,leg)),cols_gam),length(leg))
leg<-paste0(c("\u2264",rep("",length(leg)-1)),round(leg,0))
legend("bottomright",legend=rev(leg),fill=rev(col),cex=1,border=NA,bg="grey90",box.lwd=NA,inset=c(0.05,0.1),title=paste("GAM group size","\n",100*tau,"%"))
plot(land,border="grey75",add=TRUE)

#plot(g1,border="green",add=TRUE)
#plot(g2,border="red",add=TRUE)

######################################
### KERNELS KDE
######################################

#H<-Hpi.diag(coordinates(xs))
#H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
#H1<-H*matrix(c(0.02,0,0,0.02),nrow=2) 
H1<-matrix(c(6000000,0,0,6000000),nrow=2) 

### get weights from effort
ke<-kde(x=coordinates(ys),binned=TRUE,eval.points=coordinates(xs),compute.cont=FALSE,H=H1)
xs$we<-1-(ke$estimate/max(ke$estimate))

### GET POLYGONS
k<-kde(x=coordinates(xs),binned=TRUE,bgridsize=c(500,500),compute.cont=TRUE,H=H1,w=xs$nb*xs$we)
kp<-list()
perc<-c(25,50,75,95)
percw<-c("very high","high","medium","low")
trans<-c(0.9,0.7,0.5,0.3)
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
plot(land,col="grey95",border="grey75",add=TRUE)
plot(kp,add=TRUE,col=alpha("red",trans),border=NA)
mtext(paste(group,paste(month,collapse="_")),3,line=-2,font=2,adj=0.9)
#plot(xs,add=TRUE,pch=1,col=alpha("black",0.25),cex=15*nb/max(nb))
legend("bottomright",fill=alpha("red",trans),legend=paste(perc,"%"),border=NA,cex=1,bg="grey90",box.lwd=NA,inset=c(0.05,0.1),title="Kernel Contours")
#
#

### BARPLOT ABUNDANCE
par(mar=c(2,2,0,0),new=TRUE)
h<-hist(xs$nb,breaks=seq(0,max(xs$nb)+20,by=20),plot=FALSE)
h<-hist(xs$nb,breaks=seq(0,max(xs$nb)+20,by=20),col=alpha("red",0.5),border=ifelse(h$counts,"red",NA),xaxt="n",yaxt="n",)
axis(1)
axis(2)
#
#
#


####################################################################
### FINAL KERNEL OUTPUT
####################################################################

#png(paste0("D:/ebird/",paste0(group,season),"_2.png"),width=12,height=8,units="in",res=600,pointsize=14)

windows(w=21,h=12)
#par(mar=c(8.5,6,0,0),mfrow=c(1,1))
par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(xs,col="white")
shp<-list.files("D:/ebird/kernels",pattern=".shp")
layer<-paste(group,season,"ecsas",sep="_")
layer<-if(any(shp==layer)){layer}else{NULL}
if(!is.null(layer)){
  sea<-readOGR("D:/ebird/kernels",layer=layer)
  plot(gIntersection(bbox2pol(proj4string=CRS(prj)),sea,byid=TRUE),col=alpha(cols_kern,0.6),border=NA,add=TRUE)
}
plot(land,col="grey95",border="grey75",add=TRUE,lwd=0.5)

### KERNS
b<-bbox2pol(proj4string=proj4string(coastw))
#plot(gIntersection(coast,b,byid=TRUE),col=alpha("green",0.15),border=NA,add=TRUE)
#plot(kp[[i]],add=TRUE,col=cols[i],border=NA)
plot(gIntersection(coast,kp,byid=TRUE),add=TRUE,col=alpha(cols_kern,0.6),border=NA)
info<-paste0("group: ",group,"\nmonths: ",paste(month,collapse="_"),"\ndata source: EBIRD, ECSAS","\ncoast buffer: 5km")
mtext(info,side=3,line=-4,font=2,adj=0.05)
#plot(xs,add=TRUE,pch=1,col=alpha("black",0.25),cex=15*nb/max(nb))
#plot(gUnaryUnion(sea),col=alpha("blue",0.2),add=TRUE)
legend("topleft",title="Kernel Contours (risk)",fill=alpha(cols_kern,0.6),legend=paste(perc,"%","(",percw,")"),border=NA,cex=1,,bg="grey90",box.lwd=NA,inset=c(0.05,0.2))
#dev.off()


### HISTOGRAM
h<-lapply(seq_along(kp),function(i){
  o<-over(xs,kp[i,])
  res<-xs$nb[!is.na(o)]
  brks<-c(0,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000)
  brks<-brks[brks<=max(xs$nb)]
  res<-cut(res,breaks=brks,dig.lab=10)
  table(res)
})

### PLOT HISTOGRAM
par(new=TRUE,mar=c(8,6,0,0))
barplot(do.call("rbind",h),col=alpha(cols_kern,0.25),beside=TRUE,border=NA,xlab="",ylab="",las=2)
mtext("Number of records",side=2,line=4)
mtext("Abundance class",side=1,line=6.5)



### kp tot disaggregated
test<-disaggregate(gIntersection(coast,kp,byid=TRUE))
test<-do.call("rbind",test)
testp<-SpatialPoints(coordinates(test),proj4string=CRS(proj4string(test)))
o<-over(testp,g)
test<-test[!is.na(o),]
pk<-predict(fit,data.frame(lon=coordinates(test)[,1],lat=coordinates(test)[,2],stringsAsFactors=FALSE))[,1]
pk<-ifelse(pk<0,0,pk)

### GAMS predictions in KERNELS
windows(w=21,h=12)
par(mar=c(8.5,6,0,0),mfrow=c(1,1))
plot(xs,col="white")
plot(land,col="grey95",border="grey75",add=TRUE)
plot(test,border=NA,col=colo.scale(f(pk),cols_gam),add=TRUE,lwd=2)
mtext(paste(group,paste(month,collapse="_")),side=3,line=-2,font=2,adj=0.95)
leg<-finv((seq(f(1),f(max(pk)),length.out=12)))
col<-tail(colo.scale(f(c(pk,leg)),cols_gam),length(leg))
leg<-paste0(c("\u2264",rep("",length(leg)-1)),round(leg,0))
legend("bottomright",title=paste0("GAM group size at\n",paste0(100*tau,"%")," quantile predicted\n","in kernel polygons"),legend=rev(leg),fill=rev(col),cex=1,border=NA,bg="grey90",box.lwd=NA,inset=c(0.05,0.00))


### write shapefile
writeOGR(kp,dsn="D:/ebird/kernels",layer=paste(group,season,"ebird",sep="_"),driver="ESRI Shapefile",overwrite_layer=TRUE)

test<-gUnion(kp[1,],sea[1,])


### LEAFLET
kpc<-gIntersection(coast,kp,byid=TRUE)
leaflet() %>%
	addProviderTiles("Esri.WorldImagery",options = providerTileOptions(noWrap = TRUE)) %>%
	addPolygons(data=spTransform(kpc[1,],ll),stroke=FALSE,fillColor=cols[1],weight=0,fillOpacity=0.7,opacity=0) %>% 
 addPolygons(data=spTransform(kpc[2,],ll),stroke=FALSE,fillColor=cols[2],weight=0,fillOpacity=0.7,opacity=0) %>% 
 addPolygons(data=spTransform(kpc[3,],ll),stroke=FALSE,fillColor=cols[3],weight=0,fillOpacity=0.7,opacity=0) %>% 
 addPolygons(data=spTransform(kpc[4,],ll),stroke=FALSE,fillColor=cols[4],weight=0,fillOpacity=0.7,opacity=0)
	#addPolygons(data=na,stroke=FALSE,fillColor="white",weight=0,fillOpacity=0.5,opacity=0)
 #addPolylines(data=co)
	#addCircleMarkers(data=spTransform(xs[1:1000,],ll),stroke=FALSE,fillColor=cols[4],weight=0,fillOpacity=0.7,opacity=0)





















































