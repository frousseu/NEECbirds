
library(data.table)
library(RODBC)


### GET SPECIES NAMES FORM EBIRD
# needs 64 bit R
# get the e object (which is total ebird data) from the ebird.R file
g<-unique(e[,c("sp","taxo")])
g<-g[order(as.numeric(g$taxo)),]
g$db<-"ebird"
fwrite(g,"bird_groups_ebird.csv",row.names=FALSE,sep=";")


### GET SPECIES NAMES FORM ECSAS
# needs 32 bit R
ebird<-fread("bird_groups_ebird.csv",sep=";")
db<-odbcConnectAccess("D:/ebird/Master ECSAS v 3.46.mdb")
sp<-sqlFetch(db,"tblSpeciesInfo",as.is=TRUE) #?a plante et ne sait pas pourquoi
odbcClose(db)
sp<-sp[sp$Class=="Bird",]
sp$sp<-sp$English
sp$db<-"ecsas"
ecsas<-sp[,c("sp","db")]

bg<-merge(ebird[,c("sp","taxo")],ecsas[,c("sp"),drop=FALSE],all=TRUE)
bg<-bg[order(as.numeric(bg$taxo)),]

fwrite(bg,"bird_groups.csv",row.names=FALSE,sep=";")
