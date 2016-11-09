library(ff)
library(data.table)
library(sp)
library(rgdal)
library(rgeos)

#con<-file("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt")
#open(con)
#d<-read.table(con,skip=5000000,nrow=1,sep="\t")

ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=22 +datum=NAD83 +ellps=GRS80"
laea<-"+proj=laea +lat_0=50 +lon_0=-65"

keep<-c("COMMON NAME","TAXONOMIC ORDER","STATE","LATITUDE","LONGITUDE","OBSERVATION DATE","OBSERVATION COUNT","CATEGORY")


e<-fread("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",encoding = "UTF-8", na.strings = "",header = TRUE, sep = "\t", 
									colClasses = rep("character", 44),  # For simplicity
									verbose = TRUE,autostart=10,nrows=350000,select=keep)
names(e)<-gsub(" ","_",keep[match(names(e),keep)])

x<-e[!e$STATE%in%c("Alberta","Saskatchewan"),]
x<-x[CATEGORY%in%c("species"),]


xs<-SpatialPoints(matrix(as.numeric(c(x$LONGITUDE,x$LATITUDE)),ncol=2),proj4string=CRS(ll))





h<-readLines("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",n=1)
h<-strsplit(h,"\t")
h<-h[[1]]
r<-readLines("M:/SCF2016_FR/ebird/ebd_CA_relAug-2016/ebd_CA_relAug-2016.txt/ebd_CA_relAug-2016.txt",n=2)
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









