


path<-"D:/ebird/kernels"

### subset on certain shapefiles here
x<-list.files(path,pattern=".shp")
x<-x[grep("_ebird",x)]
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

plot(res[res$id=="k50",],col=alpha("red",0.15),border=NA)

res$season<-as.character(res$season)

months<-formatC("")

gsub("(.{2})", "\\1 ",res$season)




