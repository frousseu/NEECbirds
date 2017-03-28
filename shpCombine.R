library(rgdal)
library(rgeos)
library(sp)
library(FRutils)
library(ggplot2)
library(raster)
library(cleangeo)
library(data.table)
library(ecapputils)

rm(list = ls())

ROOT_DIR <- "c:/dev/NEEC"
# Directory containing the original kernel shapefiles
SHP_DIR <- file.path(ROOT_DIR, "NEECKernels")
# Output directory to produce pdf
DEST_DIR <- file.path(ROOT_DIR, "dest")
# Directory containing background maps
MAPS_DIR <- file.path(ROOT_DIR, "maps")

PROJ <-
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

ll <-
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj <-
  "+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
laea <- "+proj=laea +lat_0=50 +lon_0=-65"


###
### Clean and combine background maps
###
qc <- readOGR(MAPS_DIR, "Quebec")
atl <- readOGR(MAPS_DIR, "Atlantic")
qc <- spTransform(qc, CRS(prj))
atl <- spTransform(atl, CRS(prj))

# Clean atlatic shapefile
atl.clean <- clgeo_Clean(atl, strategy = "BUFFER", verbose = TRUE)
# Combine qc and atl
fullMap <- gUnion(qc, atl.clean)
# Save resulting shp
fm <- SpatialPolygonsDataFrame(fullMap,
                             data = data.frame(id = getSpPPolygonsIDSlots(fullMap)))
writeOGR(fm, dsn = MAPS_DIR, layer = "fullMap", driver = "ESRI Shapefile", overwrite_layer=TRUE)

## Load full map
# fullMap <- readOGR(dsn = MAPS_DIR, layer = "fullMap")


##
## Create map buffer for future use
##
buffer <- gBuffer(fullMap, width = 5000, byid = TRUE)
buffdiff <- gDifference(buffer, fullMap, byid = TRUE)
# Save resulting shp
df <- data.frame(id = getSpPPolygonsIDSlots(buffdiff))
row.names(df) <- getSpPPolygonsIDSlots(buffdiff)
buff <- SpatialPolygonsDataFrame(buffdiff,
                             data = df)
writeOGR(buff, dsn = MAPS_DIR, layer = "bufferCoast", driver = "ESRI Shapefile", overwrite_layer=TRUE)

# buffer <- readOGR(dsn = MAPS_DIR, layer = "bufferCoast")


### Load kernels 

ks <- readOGR(SHP_DIR, "NEECkernels")

# Do not keep 30% contours
kernels <- subset(ks, (!id %in% c("k20", "k30")))
kernels$uid <- paste(kernels@data$group, kernels@data$season, kernels@data$id, sep = "-")

# Combine polygons based on group, season and kernel layer
ids <- as.factor(kernels$uid)
k2 <- unionSpatialPolygons(kernels, ids)

# Groups to which we apply the buffer
bufferGroups <- c("shorebirds_waders", "waterfowl_dabbling")


# Remove observations on land
noCoastPolys <- vector(mode = "list", length = length(k2))
for (i in 1:length(k2)) {
  poly <- k2[i,]
  polId <- poly@polygons[[1]]@ID
  group <- strsplit(polId, "-")[[1]][1]
  print(sprintf("subsetting polygon: %s", polId))
  if (group %in% bufferGroups) {
    # buffer already removes coast
    print("extracting buffer")
    poly2 <- extractBuffer(poly, buffer)
  } else if (gIntersects(fullMap, poly, byid = TRUE)) {
    # if intersects with land, remove it
    print("removing land")
    poly2 <- gDifference(poly, fullMap, byid = TRUE)
    poly2@polygons[[1]]@ID <- polId
  } else {
    # do nothing
    poly2 <- poly
  } 
  noCoastPolys[[i]] <- poly2
}
noCoastKernels <- do.call(rbind, noCoastPolys)


# Update polygon ids to be able to join them to the data
nk <- noCoastKernels
for (i in seq_along(nk@polygons)) {
  newid <- gsub("\\s+.*", "", nk@polygons[[i]]@ID)
  nk@polygons[[i]]@ID <- newid
}

# Create ids in the data as well
dk <- as.data.table(kernels)
dk2 <- unique(dk, by = c("group", "season", "id"))
setkey(dk2, uid)
rownames(dk2) <- dk2$uid
newKernels <- SpatialPolygonsDataFrame(nk, dk2)


# Remove useless columns
newKernels2 <- newKernels[, -c(4,5,6, 20)]

# Write shapefile
writeOGR(newKernels2, dsn = SHP_DIR, layer = "NEECkernels_combined", driver = "ESRI Shapefile", overwrite_layer=TRUE)



# newKernels <- readOGR(dsn = SHP_DIR, layer = "kernels_combined")


# Optional : print a pdf to check the shapefile
# cols_kern <- c("red", "orange", "yellow")
# sub <- subset(newKernels, (group == "seabirds_alcids" & season == "08091011"))
# 
# printAll <- function(kernels, dest) {
#   kernels$uid <- paste(kernels$group, kernels$season, sep = "_")
#   pdf(file.path(dest, "kernels_all_layers.pdf"), width = 7, height = 7)
#   for (i in unique(kernels$uid)) {
#     sub <- subset(kernels, (uid == i))
#     
#     # par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))  
#     plot(bbox2pol(sub), border = "white")
#     land <- gIntersection(fullMap, bbox2pol(sub), byid = TRUE)
#     
#     plot(
#       land,
#       col = "grey95",
#       border = "grey75",
#       add = TRUE,
#       lwd = 0.5
#     )
#     plot(sub,
#          add = TRUE,
#          col = alpha(cols_kern, 0.6),
#          border = NA)
#     title(i)
#     
#   }
#   dev.off()
# }
# 
# 
# 
# printAll(newKernels, DEST_DIR)


