
library(OpenStreetMap)
library(sp)
library(rgdal)
library(plyr)
library(rgeos)
library(splancs)
library(dplyr)
library(tidyr)appartement
library(rgeos)
library(RColorBrewer)
library(mapplots)
library(scales)
library(plotGoogleMaps)
library(splancs)
library(tidyr)
library(dplyr)
library(svMisc)
library(tidyr)
library(raster)
library(quantreg)
library(tripack)
library(akima)
library(spatstat)
library(ks)

load("W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Code_R/shiny/urgencesapp/www/urgencesapp.RData")
d<-d[!is.na(d$Abundance),]
rm(list=ls()[ls()!="d"])
gc();gc();

path<-"W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Code_R/NEEC"

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

colo<-function(x=0:1,rescale=TRUE,col=c("blue","green","yellow","orange","red","darkred","purple"),alpha=0.5){
  if(rescale){
    x<-x/max(x)
  }
  palette<-colorRamp(col)(x)
  palette<-rgb(palette[,1],palette[,2],palette[,3],maxColorValue=256)
  alpha(palette,alpha)
}

colo.scale<-function(x,cols=c("blue","green","yellow","orange","red","purple"),center=TRUE,alpha=1){
  w<-which(is.na(x))
  if(any(w)){
    y<-x[-w]
  }else{
    y<-x
  }
  
  re<-function(a){
    if(any(w)){
      ans<-rep(NA,length(x))
      ans[which(!is.na(x))]<-a
      ans
    }else{
      a
    }
  }
  
  if(length(y)==1){
    colop<-colorRampPalette(cols)
    return(re(colop(y)))
  }  
  if(class(y)=="character"){
    colop<-colorRampPalette(cols)
    color<-colop(length(unique(y)))
    return(re(color[match(y,unique(y))]))
  }else{  
    if(all(y>=0 & y<=1)){
      color<-rgb(colorRamp(cols)(y),maxColorValue=256)
      return(re(color))
    }else{
      if(any(y<0) && center){
        m<-which.max(c(abs(min(y)),max(y)))     
        sca<-0.5/ifelse(m==1,abs(min(y)),max(y))     
        xx<-sca*y+0.5
        color<-rgb(colorRamp(cols)(xx),maxColorValue=256) 
        return(re(color))
      }else{
        color<-rgb(colorRamp(cols)((max(y)-y)/(max(y)-min(y))),maxColorValue=256) 
        return(re(color))
      }
    }
  }
}


legend.col <- function(col, lev){
  opar <- par
  n <- length(col)  
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
  bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  xx <- rep(box.cx, each = 2)
  par(xpd = TRUE)
  for(i in 1:n){
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
    box.cy[1] + (box.sy * (i)),
    box.cy[1] + (box.sy * (i)),
    box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = NA)
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
  ylim = c(min(lev), max(lev)),
  yaxt = "n", ylab = "",
  xaxt = "n", xlab = "",
  frame.plot = FALSE)
  axis(side = 4, las = 2, tick = TRUE, line = .25)
  par <- opar
}

find.nearest<-function(x,y){
  sapply(x,function(i){ 
    which.min(abs(i-y))
  })
}

mrc.shp <- readOGR(dsn = "W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Data_output/Quebec",
                   layer = "mrc_s",encoding="UTF-8")

municip.shp <- readOGR(dsn = "W:/PIU/Projet_BanqueDeDonnees/PROJET_R/Data_output/Quebec", 
                       layer= "munic_s",encoding="UTF-8") 



projected<-"+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"
proj4string(d)<-CRS(proj4string(mrc.shp))
d<-d[,c("Nom_FR","Month","Abundance","Groupe_FR")]
dd<-d
d<-spTransform(d,CRS(projected))
#dd<-spTransform(d,CRS(proj4string(fleuve)))
d<-d[!is.na(d$Month),]
municip.shp.proj<-spTransform(municip.shp,CRS(projected))
mrc.shp.proj<-spTransform(mrc.shp,CRS(projected))
fleuve<-municip.shp[is.na(municip.shp$MUS_NM_MUN),]
terre<-municip.shp[!is.na(municip.shp$MUS_NM_MUN),]
fleuve.proj<-municip.shp.proj[is.na(municip.shp.proj$MUS_NM_MUN),]
terre.proj<-municip.shp.proj[!is.na(municip.shp.proj$MUS_NM_MUN),]

windows()
par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(mrc.shp.proj)
plot(fleuve.proj,add=TRUE)
region <- getpoly()
region<-region[c(1:nrow(region),1),]
region<-SpatialPolygons(list(Polygons(list(Polygon(region)),ID=1)),proj4string=CRS(proj4string(municip.shp.proj)))
b<-bbox(region)
s<-5000
g<-GridTopology(c(b[1,1],b[2,1]),c(s,s),c(ceiling((b[1,2]-b[1,1])/s),ceiling((b[2,2]-b[2,1])/s)))
g<-SpatialGrid(g)
grid<-g
grid<-as(grid,"SpatialPolygons")
grid<-SpatialPolygonsDataFrame(grid,data=data.frame(id=1:length(grid)),match.ID=FALSE)
proj4string(grid)<-CRS(projected)
grid$id<-paste0("g",grid$id)
o<-apply(over(grid,fleuve.proj),1,function(i){all(is.na(i))})
grid<-grid[!o,]

terre.proj<-gUnaryUnion(terre.proj)
g<-gDistance(grid,terre.proj,byid=TRUE)[1,]
grid$dist<-g[match(grid$id,names(g))]



#################
### extract DAP
#################


extract_obs_dap<-function(
  
  x=d,                                # donn?es spatiales
  grid                              # polygone fourni o? trouver l'info                         
                         # taille du buffer autour de la r?gion s?lectionn?e pour cr?er le calendrier EPOQ                           
  
  ){

  o<-over(x,grid)
  x$id<-o$id
  x<-x[!is.na(x$id),]
  x<-x[!is.na(x$Abundance),]
  x<-x[!is.na(x$Month),]
  ids<-unique(x$id)
  epoq<-x[x$Base=="EPOQ",]
  cell.list<-vector(mode="list",length=length(ids))
  names(cell.list)<-ids
  for(i in seq_along(ids)){
 
   ### produce the summary
   temp<-x@data[x$id==ids[i],]
   
   #summ<-ddply(temp,.(Month,Nom_FR),function(j){
   #  a<-with(j,c(max(Abundance),round(quantile(Abundance,0.9),0),round(mean(Abundance),0),round(quantile(Abundance,0.5),0)))
   #  a
   #})
   
   summ <- temp %>% group_by(Month,Groupe_FR) %>% summarise(max=max(Abundance),q90=round(quantile(Abundance,0.9),0),mean=round(mean(Abundance),0),q50=round(quantile(Abundance,0.5),0))
  
   names(summ)[3:6]<-c("max","q90","mean","q50")
   summ$cell<-ids[i]
   
   cell.list[[i]]<-summ
   progress(i,length(ids))
 }
 ans<-do.call("rbind",cell.list)
 ans<-ans[order(ans$cell,ans$Month,-ans$max),]
 #ans<-ans[,c("cell",setdiff(names(ans),"cell"))] ### VERIFY THAT IT IS WORKING
 ans
 
}


test<-extract_obs_dap(x=d,grid=grid)

x<-test

grid.ll<-spTransform(grid,CRS(proj4string(mrc.shp)))
b<-bbox(grid.ll)
map.osm<-openmap(c(b[2,2],b[1,1]),c(b[2,1],b[1,2]),type="osm-bbike")
#map.osm2<-openproj(map.osm,proj4string(grid))
#map.osm3<-openproj(map.osm,proj4string(grid.ll))
  

#### grille
sp_grid<-function(x,grid,sp="Canard colvert",month="10",val="q90",title="",leg.title="Group size"){

  y<-as.data.frame(x[x$Month%in%month & x$Groupe_FR%in%sp,])
  col<-colo(y[,val],0.5)
  m<-match(grid$id,y$cell)
  #m<-ifelse(is.na(m),"grey85",m)
  grid$col<-col[m]
  grid$nb<-y[,val][m]
  
  #windows(width=18,height=12)
  #par(mar=c(0,0,0,0))
  plot(map.osm,raster=TRUE)
  plot(spTransform(grid,osm()),col=grid$col,border=NA,add=TRUE)
  #s<-sort(unique(y[,val]))
  #legend("right",legend=s,fill=col[match(s,y[,val])],border=NA,cex=0.75,bty="n",ncol=3)
  
  leg<-exp(seq(log(2),log(max(y[,val])),length.out=12))
  col<-tail(colo(c(y[,val],leg),0.5),length(leg))
  nb<-round(leg,0)
  #keep<-find.nearest(,nb)
  legend("left",col=rev(col),pch=15,pt.cex=2.5,legend=rev(nb),border=NA,bty="n",inset=c(0.05,0),title=leg.title)
  mtext(paste(sp,month,collapse=" "),3,line=-2,font=2)
  
  #gridk<-grid[!is.na(grid$nb),]
  #k<-kriging(coordinates(gridk)[,1],coordinates(gridk)[,2],gridk$nb,model="spherical")
  #kshp<-SpatialPoints(k$map,CRS(proj4string(grid)))
  #o<-over(kshp,grid)
  #k$map<-k$map[!is.na(o$id),]
  #windows()
  #par(mar=c(0,0,0,0))
  #plot(map.osm2)
  #image(k,col=colo(sort(unique(gridk$nb)),alpha=0.5),add=TRUE)
  grid
}


sp<-unique(d$Nom_FR[d$Groupe_FR=="Canards barboteurs"])
#sp<-"Canard noir"
month<-"05"
val<-"q90"

#o<-over(d[d$Nom_FR%in%sp & d$Month%in%month,],grid)
#table(is.na(o[,1]))

png("im0.png",width=12,height=8,units="in",res=300)
#windows(width=18,height=12)
plot(map.osm)
#plot(spTransform(grid,osm()),add=TRUE,border="grey90")
plot(spTransform(d[d$Nom_FR%in%sp & d$Month%in%month,],osm()),pch=16,col=alpha("red",0.05),add=TRUE,cex=0.6)
dev.off()


#######################################
### sp grid
#######################################

png("im1.png",width=12,height=8,units="in",res=300)
#windows(width=18,height=12)
sp_grid(x,grid,sp=sp,month=month,val=val)
dev.off()

########################################## 
### rqss
##########################################

dfit<-d[d$Nom_FR%in%sp & d$Month%in%month,]
Long<-coordinates(dfit)[,1]
Lat<-coordinates(dfit)[,2]
fit <- rqss(Abundance ~ qss(cbind(Long,Lat),lambda=0.05),data=dfit@data,tau=0.90)

g<-gBuffer(gConvexHull(dfit),width=-1000)
o<-over(SpatialPoints(coordinates(grid),CRS(proj4string(dfit))),g)
k<-!is.na(o)
X<-coordinates(grid)[k,1]
Y<-coordinates(grid)[k,2]
p<-predict(fit,data.frame(Long=X,Lat=Y,stringsAsFactors=FALSE))[,1]
p<-ifelse(p<0,0,p)

#windows(width=18,height=12)
png("im2.png",width=12,height=8,units="in",res=300)
par(mar=c(0,0,0,0))
plot(map.osm,raster=TRUE)
#plot(SpatialPoints(cbind(X,Y)),col=alpha(colo.scale(p/max(p))),pch=16,cex=3*(p/max(p)),add=TRUE)
plot(spTransform(grid[k,],osm()),col=alpha(colo.scale(p/max(p)),0.5),border=NA,add=TRUE)
#plot(spTransform(dfit,osm()),add=TRUE,pch=1)
#plot(g,add=TRUE)
leg<-exp(seq(log(2),log(max(p)),length.out=20))
col<-tail(colo(c(p,leg)/max(p),0.5),length(p))
nb<-round(leg,0)
legend("left",col=rev(col),pch=15,pt.cex=2.5,legend=rev(nb),border=NA,bty="n",inset=c(0.05,0),title="No. of individuals\nat the 90% quantile")
dev.off()


############################################
#### kernel
############################################

###kernel weighted on abundance at each point - example
#FBolduc, 01 avril 2016

d2 <- dd[which(dd$Nom_FR%in%sp),]
d2 <- d2[which(d2$Month >= month & d2$Month <= month),]
d2 <- d2[!is.na(d2$Abundance),]

#plot point pattern
hppp <- ppp( coordinates(d2)[,1],coordinates(d2)[,2],xrange=c(b[1,1],b[1,2]),yrange=c(b[2,1],b[2,2]),marks=d2@data[,"Abundance"])
spatstat.options(npixel=c(1000,1000))
birds <- density(hppp,sigma=0.03,weights=marks(hppp),positive=TRUE)

### r

png("im3.png",width=12,height=8,units="in",res=300)
m<-t(birds$v)
dat1=list()
dat1$x=seq(birds$xcol[1],by=birds$xstep,len=birds$dim[1])
dat1$y=seq(birds$yrow[1],by=birds$ystep,len=birds$dim[2])
dat1$z=m
r<-raster(dat1)  # raster complet
r<-projectRaster(r,crs=osm())
#windows(width = 18,height = 12)
plot(map.osm,raster=TRUE)
plot(r,add=TRUE,alpha=0.6,inset=c(0.2,0)) 
plot(r,legend.only=TRUE,legend.args=list(inset=c(0.2,0))) 
#plot(fleuve,add=TRUE,col="grey10")
dev.off()

### r2
png("im4.png",width=12,height=8,units="in",res=300)
m<-t(birds$v)
val<-quantile(m,0.997)
m[]<-ifelse(m<val,NA,val)
dat1=list()
dat1$x=seq(birds$xcol[1],by=birds$xstep,len=birds$dim[1])
dat1$y=seq(birds$yrow[1],by=birds$ystep,len=birds$dim[2])
dat1$z=m
r2<-raster(dat1) #raster hotspot
r2<-projectRaster(r2,crs=osm())
#windows(width = 18,height = 12)
plot(map.osm,raster=TRUE)
plot(r2,add=TRUE,alpha=0.6,col=heat.colors(100)[1],legend=FALSE) 
dev.off()
#plot(fleuve,add=TRUE)
#pr2<-rasterToPolygons(r2,dissolve=TRUE)
#plot(pr2,border="red",add=TRUE,lwd=2)

#con<-contour(birds,add=TRUE,nlevels=1)

#r3<-projectRaster(r1, crs=proj4string(d))
#r4<-projectRaster(r1, crs=osm())
#windows(width = 18,height = 12)
#plot(map.osm2)
#plot(r2,add=TRUE,alpha=0.5,col=pal_trans2)

#dev.off()


############################################
#### GROUPS
############################################

group<-c("Limicoles","Canards barboteurs")
month<-formatC(1:12,width=2,flag=0)[8:9]
case<-expand.grid(group=group,month=month,stringsAsFactors=FALSE)
val<-"q90"

############################################
#### GROUPS simple grid
############################################

png(paste0(path,"/grid.png"),width=12,height=8,units="in",res=300)
par(mar=c(0,0,0,0),mfrow=c(2,2))
l<-vector(mode="list",nrow(case))
for(i in seq_len(nrow(case))){
  l[[i]]<-sp_grid(x,grid,sp=case$group[i],month=case$month[i])
}
dev.off()

res<-lapply(1:nrow(case),function(i){
  x<-spChFIDs(l[[i]],paste(case$group[i],case$month[i],1:nrow(l[[i]]),sep="_"))
  x$group<-case$group[i]
  x$month<-case$month[i]
  x$size<-x$nb
  x
})
res<-do.call("rbind",res)
res<-res[,c("group","month","size")]

writeOGR(res[sample(1:nrow(res),100),],dsn=path,"grid",driver="ESRI Shapefile",overwrite=TRUE)

############################################
#### GROUPS rqss
############################################

png(paste0(path,"/gam_rqss.png"),width=12,height=8,units="in",res=300)
par(mar=c(0,0,0,0),mfrow=c(2,2))
l<-vector(mode="list",nrow(case))
for(i in seq_len(nrow(case))){
  sp<-unique(d$Nom_FR[d$Groupe_FR%in%case$group[i]])
  d2 <- d[which(d$Nom_FR%in%sp),]
  d2 <- d2[which(d2$Month%in%case$month[i]),]
  d2 <- d2[!is.na(d2$Abundance),]#needed for later
  w<-which(grid$dist>2000) #on met des z?ros dans les cellules ? plus de X km de la c?te
  Long<-c(coordinates(d2)[,1],coordinates(grid)[,1][w])
  Lat<-c(coordinates(d2)[,2],coordinates(grid)[,2][w])
  Abundance<-c(d2$Abundance,rep(0,length(w)))
  fit<-rqss(Abundance~qss(cbind(Long,Lat),lambda=0.05),tau=0.90)
  g<-gBuffer(gConvexHull(d2),width=-1000)
  o<-over(SpatialPoints(coordinates(grid),CRS(proj4string(d2))),g)
  k<-!is.na(o)
  X<-coordinates(grid)[k,1]
  Y<-coordinates(grid)[k,2]
  p<-predict(fit,data.frame(Long=X,Lat=Y,stringsAsFactors=FALSE))[,1]
  p<-ifelse(p<0,0,p)
  temp<-spTransform(grid[k,],osm())
  temp$group<-case$group[i]
  temp$month<-case$month[i]
  temp$size<-round(p,0)
  temp<-spChFIDs(temp,paste(temp$id,case$group[i],case$month[i],sep="_"))
  l[[i]]<-temp
  trans<-function(x,min=0.05){ans<-sqrt(x)/max(sqrt(x));ans<-ifelse(ans<min,min,ans);ans}
  plot(map.osm,raster=TRUE)
  plot(temp,col=alpha(colo.scale(p/max(p)),trans(p)),border=NA,add=TRUE)
  mtext(paste(case$group[i],case$month[i],collapse=" "),3,line=-2,font=2)
  #leg<-exp(seq(log(1),log(max(p)),length.out=12))
  leg<-(seq(sqrt(1),sqrt(max(p)),length.out=12))^2
  col<-tail(alpha(colo.scale(c(p,leg)/max(c(p,leg))),trans(c(p,leg))),length(leg))
  nb<-paste0(c("\u2264",rep("",length(leg)-1)),round(leg,0))
  legend("left",col=rev(col),pch=15,pt.cex=2.5,legend=rev(nb),border=NA,bty="n",inset=c(0.05,0),title="Group size")
  progress(i,nrow(case))
}
dev.off()

res<-do.call("rbind",l)
res<-res[,c("group","month","size")]

writeOGR(res[sample(1:nrow(res),100),],dsn=path,"gam_rqss",driver="ESRI Shapefile",overwrite=TRUE)

############################################
#### GROUPS kernel hotspots
############################################

par(mar=c(0,0,0,0),mfrow=c(1,1))
l<-vector(mode="list",nrow(case))
for(i in seq_len(nrow(case))){
  ## need to run in latlon don't know why
  sp<-unique(d$Nom_FR[d$Groupe_FR%in%case$group[i]])
  d2 <- dd[which(dd$Nom_FR%in%sp),]
  d2 <- d2[which(d2$Month%in%case$month[i]),]
  d2 <- d2[which(d2$Abundance>0),]#needed for later
  b<-bbox(fleuve)
  hppp<-ppp(coordinates(d2)[,1],coordinates(d2)[,2],xrange=c(b[1,1],b[1,2]),yrange=c(b[2,1],b[2,2]),marks=d2$Abundance)
  spatstat.options(npixel=c(1000,1000))
  h<-density(hppp,sigma=0.05,weights=marks(hppp),positive=FALSE)

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
  pr$group<-case$group[i]
  pr$month<-case$month[i]
  pr<-spChFIDs(pr,paste(case$group[i],case$month[i]))
  l[[i]]<-pr
  progress(i,nrow(case))
}



res<-do.call("rbind",l)
res<-res[,c("group","month")]
res<-spTransform(res,CRS(proj4string(d)))

png(paste0(path,"/hotspots.png"),width=12,height=8,units="in",res=300)
par(mar=c(0,0,0,0),mfrow=c(2,2))
lapply(l,function(i){
  plot(map.osm,raster=TRUE)
  mtext(paste(unlist(i@data[1,c("group","month")]),collapse=" "),3,line=-2,font=2)
  plot(spTransform(i,osm()),add=TRUE,col=alpha("red",0.5),border=NA)
})
dev.off()

writeOGR(res,dsn=path,"hotspots",driver="ESRI Shapefile",overwrite=TRUE)



###########################
### kde
###########################

d2 <- d[which(d$Nom_FR%in%sp),]
d2 <- d2[which(d2$Month >= month & d2$Month <= month),]
d2 <- d2[!is.na(d2$Abundance),]

#H<-Hpi.diag(coordinates(d2))
#H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic)
#H1<-H*matrix(c(0.5,0,0,0.5),nrow=2) 
#H2<-H*matrix(c(0.25,0,0,0.25),nrow=2) 
H2<-matrix(c(50000000,0,0,50000000),nrow=2) 
#H<-Hpi.diag(coordinates(x[x$id==ind[i],]))
k<-kde(x=coordinates(d2),compute.cont=TRUE,H=H2,w=d2$Abundance)
kp<-kde2pol(k,perc="75%",proj=proj4string(d2))

windows(width=18,height=12)
#plot(r2,xlim=unname(unlist(bbox(kp)[1,])),ylim=unname(unlist(bbox(kp)[2,])))
plot(r3)
plot(fleuve.proj,add=TRUE)
plot(kp,add=TRUE,lwd=2,border="red")
#plot(d2,pch=1,cex=(d2$Abundance/max(d2$Abundance))*5,col=alpha("blue",0.25),add=TRUE)

#o<-over(d2,kp)
#windows(width = 18,height = 12)
#par(mfrow=c(1,2))
#h1<-hist(d2$Abundance[!is.na(o)],breaks=seq(0,max(d2$Abundance)+20,by=20))
#h2<-hist(d2$Abundance[is.na(o)],breaks=seq(0,max(d2$Abundance)+20,by=20))
#barplot(rbind(h1$counts,h2$counts),names.arg=h1$mids,las=2,beside=TRUE,col=c("red","blue"),border=NA)







