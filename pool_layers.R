


path<-"D:/ebird/kernels"

x<-list.files(path,pattern=".shp")
x<-gsub("\\.shp","",x)

x<-strsplit(x,"_")
