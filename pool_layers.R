


path<-"D:/ebird/kernels"

### subset on certain shapefiles here
x<-list.files(path,pattern=".shp")
#x<-x[grep("_ebird",x)]
x<-gsub("\\.shp","",x)
x<-x[seq_along(x)]


###
lshp<-list()
for(i in seq_along(x)){
	shp<-readOGR(dsn=path,layer=x[i],encoding="UTF-8",stringsAsFactors=FALSE)
	lshp[[i]]<-shp
}

names(lshp)<-x

lshp<-lapply(seq_along(lshp),function(i){
	spChFIDs(lshp[[i]], paste(row.names(lshp[[i]]),names(lshp)[i],sep="_"))
})

res<-do.call("rbind",lshp)



res$season<-as.character(res$season)

months<-formatC(1:12,width=2,flag=0)
nmonths<-paste0("m",months)

mtable<-do.call("rbind",lapply(strsplit(gsub("(.{2})", "\\1 ",res$season)," "),function(i){as.integer(months%in%i)}))
mtable<-as.data.frame(mtable)
names(mtable)<-nmonths

slot(res,"data")<-cbind(res@data,mtable)

boxcut<-bbox2pol(c(460826,2406881,4478758,5945517),proj4string=proj4string(res))
plot(gIntersection(boxcut,res[res$id=="k50",],byid=TRUE),col=alpha("red",0.15),border=NA)


