library(data.table) #need development version 1.9.7 on github for fread to work on ebird
library(sp)
library(rgdal)
library(rgeos)
library(scales)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(spatstat)
library(raster)
library(ks)
library(FRutils)
library(svMisc)
library(RCurl)

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

### RUN THIS SCRIPT WITH 64 BIT FOR fread FUNCTION

ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
laea<-"+proj=laea +lat_0=50 +lon_0=-65"

### GET BIRD GROUPS
g<-getURL("https://raw.githubusercontent.com/frousseu/NEECbirds/master/bird_groups.csv") # Ce fichier est sur mon github
g<-read.csv(text=g,header=TRUE,stringsAsFactors=FALSE)
g<-g[!is.na(g$group_kernels),]
g<-unique(g[,c("sp","group_kernels")])
g<-g[!duplicated(g$sp),]
names(g)<-c("sp","group")

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
e$source<-"ebird"
e<-e[state%in%c("Quebec","Nova Scotia","New Brunswick","Prince Edward Island","Newfoundland and Labrador"),]

### GET SPECIES CODES
sp<-getURL("https://raw.githubusercontent.com/frousseu/ECSASatlas/master/EC_AVIAN_CORE_20161216.csv",.encoding="LATIN1") # Ce fichier est sur mon github
sp<-read.csv(text=sp,header=TRUE,stringsAsFactors=FALSE)

### ADD EPOQ DATA
load("D:/ebird/urgencesapp.RData")
rm(biomq_shp,rtlo_shp,wl,h,he,municip.shp,noga_shp,peril,razo_shp,biomq,epoqxy,perilp)
d<-d[d$Base%in%c("EPOQ"),]
l<-list(table(substr(e$date[e$state=="Quebec" & !e$group%in%c("","raptors","passerines")],1,4)),table(substr(d$Date,1,4)))
# keep all of epoq 2010 and older, based on when ebird seems to be getting more important, but there might be some overlap
d<-d[substr(d$Date,1,4)<="2010",]
d<-d@data[,c("Nom_FR","Nom_EN","Date","Month","Abundance","Long","Lat","Base","Feuillet")]
names(d)<-c("spFR","spEN","date","month","nb","lon","lat","source","sample")
m<-match(d$spFR,sp$French_Name)
d$sp<-ifelse(!is.na(m),sp$English_Name[m],NA)
d$sp<-ifelse(is.na(d$sp),d$spEN,d$sp)
d<-d[,c("spEN","date","month","nb","lon","lat","source","sample")]
names(d)<-c("sp","date","month","nb","lon","lat","source","sample")
ddply(d,.(sp),nrow) # check which species are missing
d<-d[!is.na(d$sp),] # Certains noms sp ne se retrouve pas dans les tables
d$state<-"Quebec"
d$taxo<-NA
d$date<-gsub("/","-",d$date)
d$category<-"species"
d$source<-"epoq"

### ADD BOTH DATA
x<-rbind(e,d[,names(e)])
x<-merge(x,g,all=TRUE)
x<-x[x$category%in%c("species"),] #on perd des formes et qques taxons/formes
x$month<-substr(x$date,6,7)
#e<-as.data.table(e)

#month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
#names(month_comb)<-substr(names(month_comb),1,8)

x$season<-get_season(x)
x$taxo<-as.numeric(x$taxo)
temp<-unique(x[,c("sp","taxo")])
temp<-temp[order(temp$taxo),]
x$nb<-as.numeric(x$nb)
x$nb<-ifelse(is.na(x$nb),1,x$nb)
x$lon<-as.numeric(x$lon)
x$lat<-as.numeric(x$lat)
x<-x[lat<=52,] #on garde ce qui est en bas
x<-x[lon>=(-74),] #on garde ce qui est à partir de mtl
ddply(x[x$group%in%c("",NA),],.(sp,group),nrow)




#################################################
### SEASONAL EFFORT FROM CHECKLIST LOCATIONS
#################################################

y<-unique(x[,c("sample","lat","lon","month")])
ys<-SpatialPointsDataFrame(matrix(as.numeric(c(y$lon,y$lat)),ncol=2),proj4string=CRS(ll),data=y,match.ID=FALSE)
ys<-spTransform(ys,CRS(prj))
o<-over(ys,coast)
ys<-ys[!is.na(o),]

cols_EFF<-c("blue","violet","magenta","yellow")
perc<-c(seq(5,95,by=10),99)

### GET POLYGONS
kp<-lapply(unique(get_season(x)),function(j){
	mm<-unlist(strsplit(gsub("(.{2})", "\\1 ",j)," "))
	yys<-ys[ys$month%in%mm,]	
	H1<-matrix(c(20000000,0,0,20000000),nrow=2)  
	k<-kde(x=coordinates(yys),binned=TRUE,bgridsize=c(500,500),compute.cont=TRUE,H=H1)
	kp<-list()
	trans<-rev(seq(0.15,1,length.out=length(perc)))
	col_eff<-c(colo.scale(perc[-length(perc)],cols_EFF),alpha("blue",0.5))
	kp<-kde2pol(k,levels=perc,proj4string=prj) # extract polygons
 kp$season<-j
 kp
})
names(kp)<-unique(get_season(x))


### PLOT IMAGES
# something not plotting right in the effort images
#png("D:/ebird/ebird_effort.png",width=12,height=8,units="in",res=600,pointsize=14)
par(mar=c(0,0,0,0),mfrow=c(3,2))
for(i in seq_along(kp)[1]){
  plot(bbox2pol(kp[[length(kp)]]),col="white",border="white")
  plot(land,border="grey75",add=TRUE)
  for(j in seq_along(kp[[i]])){
	   #plot(kp[[i]],add=TRUE,col=alpha("red",trans[i]),border=NA)
	   plot(gIntersection(coast,kp[[i]][j,]),add=TRUE,col=cols_EFF[j],border=NA)
  }
  legend("bottomright",fill=cols_EFF,legend=paste(perc,"%"),border=NA,cex=1,box.lwd=NA,inset=c(0.05,0.1),title="Kernel Contours")
  mtext(paste("season:",names(kp)[i]),side=3,line=-2,font=2,adj=0.9)
}
#dev.off()



#################################
### LOOP OVER GROUPS
#################################

### transform weights
f<-function(x){log(x+1)}
#f<-function(x){x}

case<-x[,.(n=.N),by=.(group,season)]
case<-case[case$n>=100 & !case$group%in%c("",NA),] # does not seem to have such a case with under 100

for(j in seq_len(nrow(case))){

  group<-case$group[j]
  #sp<-"Glaucous-winged Gull"
  #sp<-unique(x$sp[x$group%in%group])
  season<-case$season[j]
  m<-x$group%in%group & x$season%in%season
  xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll),data=x[m,])
  xs<-spTransform(xs,CRS(prj))
  o<-over(xs,coast)
  xs<-xs[!is.na(o),]


  ######################################
  ### KERNELS KDE
  ######################################

  #H<-Hpi.diag(coordinates(xs))
  #H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
  #H1<-H*matrix(c(0.02,0,0,0.02),nrow=2) 
  H1<-matrix(c(6000000,0,0,6000000),nrow=2) 

  ### get weights from effort kernel
  mm<-unlist(strsplit(gsub("(.{2})", "\\1 ",season)," "))
  ke<-kde(x=coordinates(ys[ys$month%in%mm,]),binned=TRUE,eval.points=coordinates(xs),compute.cont=FALSE,H=H1)
  xs$we<-1-(ke$estimate/max(ke$estimate))

  ### GET POLYGONS
  k<-kde(x=coordinates(xs),binned=TRUE,bgridsize=c(500,500),compute.cont=TRUE,H=H1,w=xs$nb*xs$we)
  kp<-list()
  perc<-c(25,50,75,95)
  percw<-c("very high","high","medium","low")
  trans<-c(0.9,0.7,0.5,0.3)
  cols_kern<-c("darkred","red","orange","yellow")
 	kp<-kde2pol(k,levels=perc,proj4string=proj4string(xs),cut=FALSE) # extract 
  kp$group<-group
  kp$season<-season
  kp$site<-NA
  kp$sp<-NA

  ####################################################################
  ### PLOT PNG
  ####################################################################

  png(paste0("D:/ebird/kernels/",paste(group,season,sep="_"),"_ebird.png"),width=12,height=8,units="in",res=600,pointsize=14)
  
  #par(mar=c(8.5,6,0,0),mfrow=c(1,1))
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(bbox2pol(xs),border="white")
  plot(land,col="grey95",border="grey75",add=TRUE,lwd=0.5)

  ### KERNS
  plot(gIntersection(coast,kp,byid=TRUE),add=TRUE,col=alpha(cols_kern,0.6),border=NA)
  info<-paste0("group: ",group,"\nseason: ",season,"\ndata source: EBIRD","\ncoast buffer: 5km")
  mtext(info,side=3,line=-4,font=2,adj=0.05)
  legend("topleft",title="Kernel Contours (risk)",fill=alpha(cols_kern,0.6),legend=paste(perc,"%","(",percw,")"),border=NA,cex=1,,bg="grey90",box.lwd=NA,inset=c(0.05,0.2))
#dev.off()


  ### OPTIONAL HISTOGRAM OF NUMBERS OF OBSERVATIONS
  #h<-lapply(seq_along(kp),function(i){
	 #  o<-over(xs,kp[i,])
	 #  res<-xs$nb[!is.na(o)]
	 #  brks<-c(0,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000)
	 #  brks<-brks[brks<=max(xs$nb)]
	 #  res<-cut(res,breaks=brks,dig.lab=10)
	 #  table(res)
  #})
  #par(new=TRUE,mar=c(8,6,0,0))
  #barplot(do.call("rbind",h),col=alpha(cols_kern,0.25),beside=TRUE,border=NA,xlab="",ylab="",las=2)
  #mtext("Number of records",side=2,line=4)
  #mtext("Abundance class",side=1,line=6.5)

  dev.off()
  
  ##########################################################
  ### WRITE SHAPEFILE
  ##########################################################


  ### write shapefile
  
  kp$data_type<-"ebird"
  kp$risk_level<-percw
  
  writeOGR(kp,dsn="D:/ebird/kernels",layer=paste(group,season,"ebird",sep="_"),driver="ESRI Shapefile",overwrite_layer=TRUE)


  progress(j,nrow(case))
  
}




















































