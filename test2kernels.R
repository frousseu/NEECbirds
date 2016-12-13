
H1<-matrix(c(200000000,0,0,200000000),nrow=2)

n1<-100
n2<-100

x1<-rnorm(n1,1614248,sd=20000)
y1<-rnorm(n1,5596724,sd=20000)

x2<-rnorm(n2,1409306,sd=20000)
y2<-rnorm(n2,5319450,sd=20000)

xs<-SpatialPointsDataFrame(matrix(c(x1,x2,y1,y2),ncol=2),data=data.frame(nb=c(rep(1,n1),rep(1,n2)),we=c(rep(1,n1),rep(1,n2))),proj4string=CRS(proj4string(land)))

### GET KERNELS POLYGONS
k<-kde(x=coordinates(xs),binned=TRUE,bgridsize=c(500,500),compute.cont=TRUE,H=H1,w=xs$nb*xs$we)
k2<-kde(x=coordinates(xs),binned=FALSE,eval.points=coordinates(xs),compute.cont=TRUE,H=H1,w=xs$nb*xs$we)
kp<-list()
perc<-c(25,50,75,95,99)
percw<-c("very high","high","medium","low","very low")
trans<-c(0.8,0.6,0.4,0.2,0.1)
cols_kern<-c("darkred","red","orange","yellow","lightgreen")
for(i in seq_along(perc)){
	kp[[i]]<-kde2pol(k,perc=paste0(100-perc[i],"%"),proj=proj4string(xs)) # extract polygons
}
for(i in rev(seq_along(kp))){
	if(i==1){ #make sure a single polygon for each contour
		kp[[i]]<-gUnaryUnion(kp[[i]])
	}else{	
		kp[[i]]<-gSymdifference(kp[[i]],kp[[i-1]],byid=FALSE) # keep non overlapping parts
	}
	id<-paste0("k",perc[i])
	res<-SpatialPolygonsDataFrame(kp[[i]],data=data.frame(id=id,group=group,season=season,stringsAsFactors=FALSE),match.ID=FALSE)
	kp[[i]]<-spChFIDs(res,id) # give unique ID
}
kp<-do.call("rbind",kp)

### PLOT POLYGONS
par(mar=c(1,0,0,0))
plot(xs,col="white")
plot(kp,add=TRUE,col=cols_kern,border=NA)
plot(land,col="grey95",border="grey75",add=TRUE)
legend("bottomright",fill=cols_kern,legend=paste(perc,"%"),border=NA,cex=1,bg="grey90",box.lwd=NA,inset=c(0.05,0.1),title="Kernel Contours")
plot(xs,add=TRUE,pch=1,cex=(k2$estimate/max(k2$estimate))*5)

### wr