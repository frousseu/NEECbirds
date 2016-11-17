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

na<-readOGR(dsn="D:/ebird/ebd_CA_relAug-2016.txt",layer="ne_10m_admin_1_states_provinces",encoding="UTF-8")

ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=22 +datum=NAD83 +ellps=GRS80"
laea<-"+proj=laea +lat_0=50 +lon_0=-65"

keep<-c("COMMON NAME","TAXONOMIC ORDER","STATE","LATITUDE","LONGITUDE","OBSERVATION DATE","OBSERVATION COUNT","CATEGORY","SAMPLING EVENT IDENTIFIER")
keepn<-c("sp","taxo","state","lat","lon","date","nb","category","sample")

# 29 409 415 rows in there
#e<-fread("D:/ebird/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",encoding = "UTF-8", na.strings = "",header = TRUE, sep = "\t",colClasses = rep("character", 44),verbose = TRUE,select=keep)#,nrows=29409415) 
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

x<-e[!state%in%c("Alberta","Saskatchewan"),]
x<-x[category%in%c("species"),]

x$taxo<-as.numeric(x$taxo)
temp<-unique(x[,c("sp","taxo")])
temp<-temp[order(temp$taxo),]
x$nb<-as.numeric(x$nb)
x$nb<-ifelse(is.na(x$nb),1,x$nb)
x$lon<-as.numeric(x$lon)
x$lat<-as.numeric(x$lat)
x<-x[lat<=53,] #on garde ce qui est en bas


#################################
### grouping
#################################

sp<-"Semipalmated Sandpiper"
month<-"07"
m<-x$sp%in%sp & substr(x$date,6,7)%in%month
xs<-SpatialPoints(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll))
nb<-x$nb[m]
lat<-x$lat[m]
lon<-x$lon[m]
grid<-FRutils:::hexgrid(xs,width=NULL,n=100,convex=FALSE)
#x<-x[x$taxo<=2324,]
o<-over(xs,grid)
boxplot(nb~o$id)
xx<-x[m,]




devtools:::install_github("ndphillips/yarrr")
library("yarrr")

pirateplot(nb~o$id,data=xx

)


#windows()
#par(mar=c(0,0,0,0))


#####################################
### RQSS
#####################################

fit<-rqss(nb~qss(cbind(lon,lat),lambda=0.005),tau=0.90)

g<-gBuffer(gConvexHull(xs),width=-0.001)
o<-over(SpatialPoints(coordinates(grid),CRS(proj4string(grid))),g)
w<-!is.na(o)
X<-coordinates(grid)[w,1]
Y<-coordinates(grid)[w,2]
p<-predict(fit,data.frame(lon=X,lat=Y,stringsAsFactors=FALSE))[,1]
p<-ifelse(p<0,0,p)

cols<-rev(c("white","yellow","red","darkred"))
#trans<-function(x,min=0.05){ans<-sqrt(x)/max(sqrt(x));ans<-ifelse(ans<min,min,ans);ans}
plot(xs,col="white")
plot(na,col="grey85",border="grey60",add=TRUE)
plot(grid[w,],col=colo.scale(sqrt(p),cols),border=NA,add=TRUE)
leg<-(seq(sqrt(1),sqrt(max(p)),length.out=12))^2
col<-tail(colo.scale(sqrt(c(p,leg)),cols),length(leg))
leg<-paste0(c("\u2264",rep("",length(leg)-1)),round(leg,0))
legend("topright",legend=rev(leg),fill=rev(col),cex=2)
#plot(xs,add=TRUE,pch=1,col=alpha("black",0.005))
#plot(na,add=TRUE)


######################################
### KERNELS KDE
######################################

H<-Hpi.diag(coordinates(xs))
H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic)
#H1<-H*matrix(c(0.5,0,0,0.5),nrow=2) 
H1<-H*matrix(c(0.25,0,0,0.25),nrow=2) 
#H1<-matrix(c(50000000,0,0,50000000),nrow=2) 

k<-kde(x=coordinates(xs),compute.cont=TRUE,H=H1,w=nb)
kp<-kde2pol(k,perc="75%",proj=proj4string(grid))
plot(kp,add=TRUE,col=alpha("darkgreen",0.25),border="darkgreen")



#
#
#
#
#
#
#







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
















































