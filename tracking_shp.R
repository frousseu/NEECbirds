library(data.table)
library(sp)
library(rgdal)
library(rgeos)
library(scales)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(raster)
library(ks)
library(FRutils)
library(svMisc)
library(RODBC)
library(readxl)
library(plyr)
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

d<-list()
month_comb<-unlist(list("12010203"=c("12","01","02","03"),"04050607"=c("04","05","06","07"),"08091011"=c("08","09","10","11")))
names(month_comb)<-substr(names(month_comb),1,8)

### GET BIRD GROUPS
g<-getURL("https://raw.githubusercontent.com/frousseu/NEECbirds/master/bird_groups.csv") # Ce fichier est sur mon github
g<-read.csv(text=g,header=TRUE,stringsAsFactors=FALSE)
g<-g[!is.na(g$group_kernels),]
g<-unique(g[,c("sp","group_kernels")])

### GET SPECIES CODES
sp<-getURL("https://raw.githubusercontent.com/frousseu/ECSASatlas/master/EC_AVIAN_CORE_20161216.csv",.encoding="LATIN1") # Ce fichier est sur mon github
sp<-read.csv(text=sp,header=TRUE,stringsAsFactors=FALSE)

######################################
### NOGA
######################################

x<-list.files("D:/Télémétrie/Donnees_GPS")
x<-x[-grep("ali",x)]
l<-sapply(x,function(i){
	path<-paste0("D:/Télémétrie/Donnees_GPS/",i) 
	sh<-excel_sheets(path)
})
l<-data.frame(f=rep(names(l),sapply(l,length)),s=unlist(l))
l<-apply(l,1,function(i){
	ans<-read_excel(path=paste0("D:/Télémétrie/Donnees_GPS/",i[1]),sheet=i[2])
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
x<-x[,c("id","Date","Latitude","Longitude")]
names(x)<-c("id","date","lat","lon")
x$sp<-"NOGA"
n<-ddply(x,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$group<-"seabirds_pelagics"
#x$season<-names(month_comb)[match(x$month,month_comb)]
x$season<-get_season(x)
x$site<-"Forillon"
d[["NOGA_Forillon"]]<-x


##############################################
### RTLO
##############################################

x<-read.csv("D:/Télémétrie/RTLO_EEZ.txt",stringsAsFactors=FALSE)
x$x<-x$longitud
x$y<-x$latitude
x$id<-x$animal
x$id2<-1
x<-x[,c("animal","date","latitude","longitud")]
names(x)<-c("id","date","lat","lon")
x$date<-substr(strptime(x$date,"%m/%d/%Y %H:%M:%S"),1,10)
x$sp<-"RTLO"
n<-ddply(x,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$group<-"waterfowl_diving"
#x$season<-names(month_comb)[match(x$month,month_comb)]
x$season<-get_season(x)
x$site<-"Unknown"
d[["RTLO_Unknown"]]<-x



###########################################
### RAZO
###########################################

#db<-odbcConnectAccess2007("D:/Télémétrie/ID_General_RAZO.accdb")
#d2015<-sqlFetch(db,"Donnees_2015",stringsAsFactors=FALSE)
#d2016<-sqlFetch(db,"Donnees_2016",stringsAsFactors=FALSE)
#odbcClose(db)
d2015<-as.data.frame(fread("D:/Télémétrie/Donnees_2015_RAZO.txt",stringsAsFactors=FALSE,na.strings=c("","NA")))
d2016<-as.data.frame(fread("D:/Télémétrie/Donnees_2016_RAZO.txt",stringsAsFactors=FALSE,na.strings=c("","NA")))
names(d2015)[which(names(d2015)=="Longtitude")]<-"Longitude"
x<-full_join(d2015,d2016)
x<-x[!is.na(x$Latitude),]
x$Longitude<-as.numeric(x$Longitude)
x$Latitude<-as.numeric(x$Latitude)
x$x<-x$Longitude
x$y<-x$Latitude
x$id<-x$LoggerID
x$id2<-1
x<-x[x$id%in%names(table(x$id)[table(x$id)>2]),]
x<-x[,c("LoggerID","date","Latitude","Longitude")]
names(x)<-c("id","date","lat","lon")
x$date<-substr(strptime(x$date,"%m/%d/%Y %H:%M:%S"),1,10)
x$sp<-"RAZO"
n<-ddply(x,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$group<-"seabirds_alcids"
#x$season<-names(month_comb)[match(x$month,month_comb)]
x$season<-get_season(x)
x$site<-"Pèlerin"
d[["RAZO_Pèlerin"]]<-x


#############################################
### GSGO
#############################################

x<-as.data.frame(fread("D:/Télémétrie/data2012_01_17_OIES.txt",stringsAsFactors=FALSE))
x<-x[x$typ=="gps" & x$xx%in%0:5 & x$good,]
x<-x[x$lat>45 & x$lat<50,]
x<-x[x$lon>(-77),]
x$x<-x$lon
x$y<-x$lat
x<-x[,c("band","date","lat","lon")]
names(x)<-c("id","date","lat","lon")
x$sp<-"GSGO"
n<-ddply(x,.(id),nrow)
n$val<-1/(n$V1/max(n$V1))
x$we<-n$val[match(x$id,n$id)]
x$month<-substr(x$date,6,7)
x$group<-"waterfowl_dabbling"
#x$season<-names(month_comb)[match(x$month,month_comb)]
x$season<-get_season(x)
x$site<-"Québec"
d[["GSGO_Québec"]]<-x


#############################################
### Atlantic tracking
#############################################

x<-as.data.frame(fread("D:/Télémétrie/TRACKING_atlantic_20170228.csv",stringsAsFactors=FALSE))
x<-x[,c("bird.id","date","lat","lon","species","site")]
names(x)<-c("id","date","lat","lon","sp","site")
x$sp<-toupper(x$sp)
x$month<-substr(x$date,6,7)
x$group<-g$group_kernels[match(sp$English_Name[match(x$sp,sp$Species_ID)],g$sp)]
#x$season<-names(month_comb)[match(x$month,month_comb)]
x$season<-get_season(x)
#x$sp<-ifelse(x$sp=="RAZO","RAZOAT",x$sp)
x2<-dlply(x,.(sp,site),function(i){
	n<-ddply(i,.(id),nrow)
	n$val<-1/(n$V1/max(n$V1))
	i$we<-n$val[match(i$id,n$id)]
	i<-i[,names(d[[1]])]
	i
})
names(x2)<-gsub("\\.","_",names(x2))
for(i in seq_along(x2)){
	d[[names(x2)[i]]]<-x2[[i]]
}

ddply(x,.(sp,site),function(i){length(unique(i$id))})
ddply(x,.(sp,site,id),nrow)
sapply(d,function(i){length(unique(i$id))})

#####################################
### BUILD SEASONAL SUBSETS
#####################################

season<-unique(x$season)
d<-unlist(lapply(d,function(i){
     res<-lapply(season,function(j){
     	 mm<-unlist(strsplit(gsub("(.{2})", "\\1 ",j)," "))
    	  ans<-i[i$month%in%mm,]
    	  if(nrow(ans)==0){
    	    "empty"
    	  }else{
    	    ans
    	  }
     })
     names(res)<-season
     res
}),recursive=FALSE)
names(d)<-gsub("\\.","_",names(d))
d<-d[sapply(d,function(i){!identical(i,"empty")})]

count<-ldply(d,function(x){
	data.frame(nb_ind=length(unique(x$id)),nb_loc=nrow(x))
})

# take out any population with fewer than 5 ids and 50 locations
keep<-count$.id[count$nb_ind>=5 & count$nb_loc>=50]
d<-d[keep]

# take out combinations that appear in 0607 and 04050607 and keep tha good one
temp<-do.call("rbind",d)
temp$season<-get_season(temp)
comb1<-apply(unique(temp[,c("sp","season")]),1,function(i){paste(i[1],i[2],sep="_")})
comb2<-sapply(strsplit(names(d),"_"),function(i){paste(i[1],i[3],sep="_")})
keep<-comb2%in%comb1
d<-d[keep]

#####################################
### BUILD SPECIES Subsets
#####################################

d2<-list()
spec<-sapply(strsplit(names(d),"_"),"[",1)
seas<-sapply(strsplit(names(d),"_"),"[",3)
case<-unique(data.frame(spec=spec,seas=seas))
for(j in seq_len(nrow(case))){
  w<-which(spec==case$spec[j] & seas==case$seas[j]) 
  d2[[paste(case$spec[j],case$seas[j],sep="_")]]<-do.call("rbind",d[w])
}
d<-d2

d<-lapply(d,function(x){
	 x$x<-x$lon
	 x$y<-x$lat
	 coordinates(x)<-~x+y
  proj4string(x)<-CRS("+proj=longlat +datum=WGS84")
  x<-spTransform(x,CRS("+proj=utm +zone=18 +datum=NAD83 +ellps=GRS80"))
  x
})

ldply(d,function(x){
	data.frame(nb_ind=length(unique(x$id)),nb_loc=nrow(x))
})


####################################################
### GET KERNELS
####################################################


for(j in seq_along(d)){
	 x<-d[[j]]
	 group<-strsplit(names(d)[j],"_")[[1]][1]
	 season<-strsplit(names(d)[j],"_")[[1]][2]
	 #H<-Hpi.diag(coordinates(x))
	 #H[H>0]<-min(H[H>0]) #this thing assumes the same variability in both directions (isotropic) and imposes the smallest value
	 #H1<-H*matrix(c(0.02,0,0,0.02),nrow=2) 
	 H1<-matrix(c(20000000,0,0,20000000),nrow=2) 
	 perc<-c(25,50,75,95)
	 percw<-c("very high","high","medium","low")
	 trans<-c(0.9,0.7,0.5,0.3)
	 cols_kern<-c("darkred","red","orange","yellow")
	 site<-unique(x$site)
	 l<-lapply(site,function(i){
	 	 xx<-x[x$site==i,]
    k<-kde(x=coordinates(xx),binned=TRUE,bgridsize=c(500,500),compute.cont=TRUE,H=H1,w=xx$we)
    kp<-kde2pol(k,levels=perc,proj4string=prj,cut=FALSE) # extract polygons
    kp<-spChFIDs(kp,paste(row.names(kp),i,sep="_"))
    kp
	 })
  kp<-do.call("rbind",l)  
  kp$group<-group
  kp$season<-season
  kp$site<-sapply(strsplit(row.names(kp),"_"),"[",2)
  
  
  ####################################################################
  ### PLOT PNG
  ####################################################################
  
  png(paste0("D:/ebird/kernels/",paste(group,season,sep="_"),"_tracking.png"),width=12,height=8,units="in",res=600,pointsize=14)
  
  #par(mar=c(8.5,6,0,0),mfrow=c(1,1))
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  if(gArea(bbox2pol(kp))>=(100000*100000)){
  	plot(bbox2pol(kp),border="white") 
  }else{
  	b<-bbox2pol(as.vector(bbox(kp))[c(1,3,2,4)]+c(-100000,100000,-100000,100000))
  	plot(b,border="white")
  }
  plot(land,col="grey95",border="grey75",lwd=0.5,add=TRUE)
  ### KERNS
  plot(kp,add=TRUE,col=cols_kern,border=NA)
  #plot(gIntersection(coast,kp,byid=TRUE),add=TRUE,col=alpha(cols_kern,0.6),border=NA)
  info<-paste0("group: ",group,"\nseason: ",season,"\ndata source: tracking")
  mtext(info,side=3,line=-4,font=2,adj=0.05)
  legend("topleft",title="Kernel Contours (risk)",fill=cols_kern,legend=paste(perc,"%","(",percw,")"),border=NA,cex=1,bg="grey90",box.lwd=NA,inset=c(0.05,0.2))
  
  dev.off()
  
  ##########################################################
  ### WRITE SHAPEFILE
  ##########################################################
  
  kp$data_type<-"tracking"
  kp$risk_level<-percw
  
  writeOGR(kp,dsn="D:/ebird/kernels",layer=paste(group,season,"tracking",sep="_"),driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  progress(j,length(d))
  
  
}



















