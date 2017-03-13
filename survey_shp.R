
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

### GET ECSAS Alpha codes

load("D:/ebird/odbcECSAS.RData") #load session instead of running odbcConnect on 32bit all the time
alpha<-sp
rm(list=ls()[ls()!="alpha"])

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
d<-d[d$Base%in%c("CANARDS_MER","EIDERS_HIVER","GARROTS_HIVER","LIMICOLES_IDLM","MACREUSES","OIES","REKN_BBPL_YT","SRIV","SAUVFLEUVE"),]
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
x<-d[d$sp%in%"Common Eider" & d$source%in%"EIDERS_HIVER",] ### tests what's left from eiders
x<-d[d$Nom_FR%in%"Eider à duvet" & d$Base%in%"EIDERS_HIVER",] ### tests what's left from eiders


### ATLANTIC DATA (from Rob Ronconi)
s<-as.data.frame(read_excel(path="D:/ebird/Sightings_data_20170228.xlsx",sheet=1),stringsAsFactors=FALSE)
source.id<-read_excel(path="D:/ebird/Sightings_data_20170228.xlsx",sheet=2)
#s<-s[apply(s,1,function(i){!all(is.na(i))}),] # des lignes vides sont présentes à la fin dans le vieux fichier excel de base
s<-s[,c("species","date","month","lat","lon","count","source.id")]
names(s)<-c("sp","date","month","lat","lon","nb","source")
s$month<-formatC(s$month,width=2,flag=0) # some missing dates so use month instead
m<-match(s$sp,sp$Species_ID)
s$sp<-ifelse(is.na(m),s$sp,sp$English[m])
#s$sp[grep("PIPL",s$sp)]<-"Piping Plover"
s$sp[grep("HARD",s$sp)]<-"Harlequin Duck"
s$sp[grep("UNPH",s$sp)]<-"Phalarope"
s$region<-"AT"

### WINTER EIDERS
# que fait-on avec les inconnus?
# on additionne les bruns et blancs
e1<-as.data.frame(read_excel(path="D:/ebird/COEIObsNBNS.xlsx",sheet=1),stringsAsFactors=FALSE)
e2<-as.data.frame(read_excel(path="D:/ebird/COEIObsQCNL.xlsx",sheet=1),stringsAsFactors=FALSE)
e1$region<-"AT"
e2$region<-"QC"
keep<-c("Mois","visuelblancs","visuelbruns","LatDec","LongDec","region","An","Mois","Jour","region")
e<-rbind(e1[,keep],e2[,keep])
e$nb<-e$visuelblancs+e$visuelbruns
#check wrong lats
e <- e[order(e$LatDec),]
#some lat-long reversed - wrong column
e$lat <- e$LatDec
e$lon <- e$LongDec
e$LongDec <- ifelse(e$LatDec<0,e$LatDec,e$LongDec)
e$LatDec <- ifelse(e$LatDec<0,e$lon,e$LatDec)
e<-e[(order(e$LatDec)),]
e$lat <- e$LatDec
e$lon <- e$LongDec
e$month<-formatC(e$Mois,width=2,flag=0)
e$date<-paste(e$An,e$mois,formatC(e$Jour,width=2,flag=0),sep="-")
e$sp<-"Common Eider"
e$source<-"Winter_Eiders"
e<-e[,c("sp","date","month","lat","lon","nb","source","region")]

### BUILD SEASONS
#month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
#names(month_comb)<-substr(names(month_comb),1,8)



### JOIN BOTH
x<-join(d,s,type="full")
x<-join(x,e,type="full")

m<-match(x$sp,g$sp)
x$sp<-ifelse(is.na(m),x$sp,g$sp[m])
m<-match(tolower(x$sp),g$sp)
x$sp<-ifelse(is.na(m),x$sp,g$sp[m])
m<-match(gsub(" sp"," sp.",tolower(x$sp)),g$sp)
x$sp<-ifelse(is.na(m),x$sp,g$sp[m])
w<-which(x$group%in%c("",NA))
x$sp[w]<-alpha$English[match(x$sp[w],alpha$Alpha)]
x$sp<-gsub(" (Larus, Xema, Rissa, Pagophila, Rhodostethia)","",x$sp)
#table(x$sp[is.na(match(x$sp,g$sp))])
x$group<-g$group[match(x$sp,g$sp)]
x$group[which(x$sp%in%c("Phalarope"))]<-"shorebirds_waders"
x$group[which(x$sp%in%c("ALCI","MURA","UNMU"))]<-"seabirds_alcids"
x$group[which(x$sp%in%c("UNSH","UNSK"))]<-"seabirds_pelagics"
x$group[which(x$sp%in%c("UNLA","UNGU","UNTE"))]<-"seabirds_larids"
x$group[which(x$sp%in%c("UNME","UNLO","UNSC"))]<-"waterfowl_diving"
table(x$sp[x$group%in%c("",NA)])
x$season<-get_season(x)

### take out what is in FB eliminate file
w1<-which(x$group%in%c("passerines","raptors") & x$source%in%c("GARROTS_HIVER","MACREUSES","SAUVFLEUVE","SRIV","CANARDS_MER"))
w2<-which(x$group%in%c("seabirds_alcids") & x$source%in%c("MACREUSES","SRIV"))
w3<-which(x$group%in%c("seabirds_alcids") & x$source%in%c("GARROTS_HIVER","CANARDS_MER") & !x$season%in%c("12010203"))
w4<-which(x$group%in%c("seabirds_larids") & x$source%in%c("GARROTS_HIVER"))
w5<-which(x$group%in%c("seabirds_pelagics") & x$source%in%c("MACREUSES","SRIV","CANARDS_MER"))
w6<-which(x$group%in%c("shorebirds_waders") & x$source%in%c("CANARDS_MER"))
x<-x[-c(w1,w2,w3,w4,w5,w6),]

### scrap zeros or missing numbers
x<-x[!is.na(x$nb) & x$nb>0,]
ddply(x[x$group%in%c("",NA),],.(sp,group),nrow)
ddply(x[x$sp%in%"Common Eider",],.(sp,season,month,region),nrow)
# check if alcids in CANARDS_MER and GARROTS_HIVER are to include
wa<-x[x$group=="seabirds_alcids" & x$source%in%c("CANARDS_MER","GARROTS_HIVER") & x$season=="12010203",]
ddply(x,.(group,source),function(i){sum(i$nb)})


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
	perc<-c(20,50,70,90)
	percw<-c("very high","high","medium","low")
	trans<-c(0.8,0.6,0.4,0.2)
	cols_kern<-c("darkred","red","orange","yellow")
	kp<-kde2pol(k,levels=perc,proj4string=proj4string(xs),cut=FALSE) # extract polygons
	kp$group<-group
	kp$season<-season
	kp$site<-NA
	kp$sp<-NA
	
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








