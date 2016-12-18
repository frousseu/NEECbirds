
library(readxl)

load("D:/ebird/urgencesapp.RData")
rm(biomq_shp,rtlo_shp,wl,h,he,municip.shp,noga_shp,peril,razo_shp,biomq,epoqxy,perilp)

d<-d[d$Base%in%c("CANARDS_MER","EIDERS_HIVER","GARROTS_HIVER","LIMICOLES_IDLM","MACREUSES","OIES","REKN_BBPL_YT","SAUVAGINE","SRIV","SAUVFLEUVE"),]



ans<-read_excel(path="D:/ebird/Sightings_data_20161213.xlsx")

