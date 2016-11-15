library(ff)
library(data.table)
library(sp)
library(rgdal)
library(rgeos)
library(quantreg)
library(FRutils)
library(scales)

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
#e<-fread("D:/ebird/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",encoding = "UTF-8", na.strings = "",header = TRUE, sep = "\t",colClasses = rep("character", 44),verbose = TRUE,select=keep)#,nrows=30000000) 

#fwrite(e[e$STATE=="Quebec",],"D:/ebird/ebd_CA_relAug-2016.txt/ebirdQC.csv",sep="\t")

e<-fread("D:/ebird/ebd_CA_relAug-2016.txt/ebirdQC.csv",encoding = "UTF-8", na.strings = "",header = TRUE, sep = "\t",verbose = TRUE)

names(e)<-keepn[match(names(e),keep)]

x<-e[!e$state%in%c("Alberta","Saskatchewan"),]
x<-x[category%in%c("species"),]
x$taxo<-as.numeric(x$taxo)
temp<-unique(x[,c("sp","taxo")])
temp<-temp[order(temp$taxo),]
x<-x[x$taxo<=2324,]
x$nb<-as.numeric(x$nb)
x$nb<-ifelse(is.na(x$nb),1,x$nb)
x$lon<-as.numeric(x$lon)
x$lat<-as.numeric(x$lat)


sp<-"Northern Shoveler"
month<-"04"
m<-x$sp%in%sp & substr(x$date,6,7)%in%month
xs<-SpatialPoints(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll))
nb<-x$nb[m]
lat<-x$lat[m]
lon<-x$lon[m]
grid<-hexgrid(xs,width=20000,convex=FALSE)

plot(xs,col="white")
#plot(grid,add=TRUE,col=alpha("red",0.25))
#g1<-!apply(gContains(gUnaryUnion(na),grid,byid=TRUE),1,any)
#grid<-grid[g1,]
fit<-rqss(nb~qss(cbind(lon,lat),lambda=0.005),tau=0.90)
#plot(fit,col="red",asp=1.5)

g<-gBuffer(gConvexHull(xs),width=-0.001)
o<-over(SpatialPoints(coordinates(grid),CRS(proj4string(grid))),g)
k<-!is.na(o)
X<-coordinates(grid)[k,1]
Y<-coordinates(grid)[k,2]

p<-predict(fit,data.frame(lon=X,lat=Y,stringsAsFactors=FALSE))[,1]
p<-ifelse(p<0,0,p)
plot(grid[k,],col=colo.scale(p,c("red","white")),border="white",add=TRUE)
plot(xs,add=TRUE,pch=1,col=alpha("black",0.005))
plot(na,add=TRUE)


xfit<-x[x$sp%in%sp,]
lon<-xfit$lon
lat<-xfit$lat
fit <- rqss(nb ~ qss(cbind(lon,lat),lambda=0.05),data=xfit,tau=0.90)

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












h<-readLines("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",n=1)
h<-strsplit(h,"\t")
h<-h[[1]]
r<-readLines("C:/Users/User/Documents/SCF2016_FR/ebird/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",n=2)
#writeLines(r,"M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/test.txt")


r<-strsplit(r,"\t")
r[[1]]<-gsub(" ","_",r[[1]])

x<-as.data.frame(do.call("rbind",r[-1]),stringsAsFactors=FALSE)
names(x)<-r[[1]]

head(x)

#scan("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",what=list(NULL),sep='\n',skip=1000000,nlines=20)

x<-read.csv.ffdf(file="M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",skip=100000,sep="\t",first.rows=30000,next.rows=30000,nrows=100000,VERBOSE=TRUE)



x<-x[x[,"CATEGORY"]=="species",c("COMMON.NAME","TAXONOMIC.ORDER","STATE","LATITUDE","LONGITUDE","OBSERVATION.DATE","OBSERVATION.COUNT")]
x<-as.data.frame(x)











chunksize=1000

con <- file("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt", "r", blocking = FALSE) #create file connection
d=scan(con,what="character",nlines=1,sep="\t") #remove the header line
se<-seq(1,3000,chunksize)
l<-vector(mode="list",length=length(se))
n<-1
for(i in se){
	
	d=scan(con,what="character",nlines=chunksize,sep="\t",quiet=TRUE)
	d=t(matrix(d,nrow=44))
	d=data.frame(d)
	l[[n]]<-d
	n<-n+1
	
	
	
	#Do stuff with d....
	
}









