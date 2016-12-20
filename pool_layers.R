


path<-"D:/ebird/kernels"

x<-list.files(path,pattern=".shp")
x<-gsub("\\.shp","",x)
x<-x[1:18]



lshp<-list()
for(i in seq_along(x)){
	shp<-readOGR(dsn=path,layer=x[i],encoding="UTF-8")
	lshp[[i]]<-shp
}

names(lshp)<-x

lshp<-lapply(seq_along(lshp),function(i){
	spChFIDs(lshp[[i]], paste(row.names(lshp[[i]]),names(lshp)[i],sep="_"))
})


res<-do.call("rbind",lshp)

plot(res[res$id=="k50",],col=alpha("red",0.25))
