
library(RCurl)
library(readxl)
library(plyr)
library(data.table)
library(ks)
library(FRutils)
library(scales)
library(sp)
library(rgeos)
library(maptools)
library(rgdal)
library(svMisc)

### DATA FROM QC AND ATL ARE SEPARATED BECAUSE QC HAS 4 TIMES MORE DATA

ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
laea<-"+proj=laea +lat_0=50 +lon_0=-65"

### GET BIRD GROUPS
g<-getURL("https://raw.githubusercontent.com/frousseu/NEECbirds/master/bird_groups.csv") # Ce fichier est sur mon github
g<-read.csv(text=g,header=TRUE,stringsAsFactors=FALSE)
g<-g[!is.na(g$group_kernels),]
g<-unique(g[,c("sp","group_kernels")])

### GET SPECIES CODES
sp<-getURL("https://raw.githubusercontent.com/frousseu/ECSASatlas/master/EC_AVIAN_CORE_20161216.csv",.encoding="LATIN1") # Ce fichier est sur mon github
sp<-read.csv(text=sp,header=TRUE,stringsAsFactors=FALSE)


### QUEBEC DATA (using the data in the UrgencesAviR workspace)
load("D:/ebird/urgencesapp.RData")
rm(biomq_shp,rtlo_shp,wl,h,he,municip.shp,noga_shp,peril,razo_shp,biomq,epoqxy,perilp)
d<-d[d$Base%in%c("CANARDS_MER","EIDERS_HIVER","GARROTS_HIVER","LIMICOLES_IDLM","MACREUSES","OIES","REKN_BBPL_YT","SAUVAGINE","SRIV","SAUVFLEUVE"),]
d<-d@data[,c("Nom_FR","Nom_EN","Date","Month","Abundance","Long","Lat","Base")]
names(d)<-c("spFR","spEN","date","month","nb","lon","lat","source")
m<-match(d$spFR,sp$French_Name)
d$sp<-ifelse(!is.na(m),sp$English_Name[m],NA)
d$sp<-ifelse(is.na(d$sp),d$spEN,d$sp)
d<-d[,c("spEN","date","month","nb","lon","lat","source")]
names(d)<-c("sp","date","month","nb","lon","lat","source")
ddply(d,.(sp,spFR,spEN),nrow) # check which species are missing
d<-d[!is.na(d$sp),] # Certains noms sp ne se retrouve pas dans les tables
d$month<-ifelse(d$month=="00","09",d$month) #there is a bug in the data for SAUVFLEUVE dealing with months and the original excel file not read correctly by read_excel
d$region<-"QC"

### ATLANTIC DATA (from Rob Ronconi)
s<-as.data.frame(read_excel(path="D:/ebird/Sightings_data_20161213.xlsx",sheet=1),stringsAsFactors=FALSE)
source.id<-read_excel(path="D:/ebird/Sightings_data_20161213.xlsx",sheet=2)
s<-s[apply(s,1,function(i){!all(is.na(i))}),] # des lignes vides sont présentes à la fin dans le fichier excel de base
s<-s[,c("species","date","month","lat","lon","count","source.id")]
names(s)<-c("sp","date","month","lat","lon","nb","source")
s$month<-formatC(s$month,width=2,flag=0) # some missing dates so use month instead
m<-match(s$sp,sp$Species_ID)
s$sp<-ifelse(is.na(m),s$sp,sp$English[m])
#s$sp[grep("PIPL",s$sp)]<-"Piping Plover"
s$sp[grep("HARD",s$sp)]<-"Harlequin Duck"
s$sp[grep("UNPH",s$sp)]<-"Phalarope"
s$region<-"AT"


### BUIDL SEASONS
month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
names(month_comb)<-substr(names(month_comb),1,8)



### JOIN BOTH
x<-join(d,s,type="full")
x$season<-names(month_comb)[match(x$month,month_comb)]

m<-match(x$sp,g$sp)
x$sp<-ifelse(is.na(m),x$sp,g$sp[m])
m<-match(tolower(x$sp),g$sp)
x$sp<-ifelse(is.na(m),x$sp,g$sp[m])
m<-match(gsub(" sp"," sp.",tolower(x$sp)),g$sp)
x$sp<-ifelse(is.na(m),x$sp,g$sp[m])
x$group<-g$group[match(x$sp,g$sp)]

### scrap zeros or missing numbers
x<-x[!is.na(x$nb) & x$nb>0,]


#################################
### grouping
#################################

### transform weights
f<-function(x){log(x+1)}
#f<-function(x){x}

case<-ddply(x,.(group,season,region),nrow)
case<-case[!case$group%in%c("passerines","raptors",NA,"") & case$V1>=100,]

for(j in seq_len(nrow(case))){
	
	group<-case$group[j]
	season<-case$season[j]
	region<-case$region[j]
	#sp<-"Dovekie"
	m<-x$group%in%group & x$season%in%season & x$region%in%region
	xs<-SpatialPointsDataFrame(matrix(as.numeric(c(x$lon[m],x$lat[m])),ncol=2),proj4string=CRS(ll),data=x[m,])
	xs<-spTransform(xs,CRS(prj))
	
	
	######################################
	### KERNELS KDE
	######################################
	
	#H<-Hpi.diag(coordinates(xs))
	#H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
	#H1<-H*matrix(c(0.02,0,0,0.02),nrow=2) 
	H1<-matrix(c(50000000,0,0,50000000),nrow=2) 
	
	### GET KERNELS POLYGONS
	k<-kde(x=coordinates(xs),binned=TRUE,gridsize=c(500,500),compute.cont=TRUE,H=H1,w=f(xs$nb))
	kp<-list()
	perc<-c(25,50,75,95)
	percw<-c("very high","high","medium","low")
	trans<-c(0.8,0.6,0.4,0.2)
	cols_kern<-c("darkred","red","orange","yellow")
	kp<-kde2pol(k,levels=perc,proj4string=proj4string(xs),cut=FALSE) # extract polygons
	kp$group<-group
	kp$season<-season
	
	png(paste0("D:/ebird/kernels/",paste(group,season,sep="_"),paste0("_survey",region),".png"),width=12,height=8,units="in",res=500,pointsize=14)
	
	### PLOT POLYGONS
	par(mar=c(1,0,0,0))
	plot(xs,col="white")
	plot(kp,add=TRUE,col=alpha(cols_kern,0.6),border=NA)
	plot(land,col="grey95",border="grey75",add=TRUE)
	legend("bottomright",fill=alpha(cols_kern,0.6),legend=paste(perc,"%"),border=NA,cex=1,bg="grey90",box.lwd=NA,inset=c(0.05,0.1),title="Kernel Contours")
	info<-paste0("group: ",group,"\nseason: ",season,paste0("\ndata source: SURVEYS",region),"\ncoast buffer:")
	mtext(info,side=3,line=-4,font=2,adj=0.05)
	
	dev.off()
	
	kp$data_type<-"survey"
	kp$risk_level<-percw
	
	### write shapefile
	writeOGR(kp,dsn="D:/ebird/kernels",layer=paste(group,season,paste0("survey",region),sep="_"),driver="ESRI Shapefile",overwrite_layer=TRUE)
	
	progress(j,nrow(case))
	
}








