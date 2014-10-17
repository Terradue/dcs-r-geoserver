#!/usr/bin/Rscript --vanilla --slave --quiet

# load rciop library to access the developer cloud sandbox functions
library("rciop")
library("rgeos")
library("stringr")

# load the application package when mvn installed it
library(rGeoServer, lib.loc="/application/share/R/library")
library(rLandsat8, lib.loc="/application/share/R/library")
load("/application/.geoserver.authn.RData")

# get the GeoServer REST access point
geoserver <- rciop.getparam("geoserver")
# get the extent of the area of interest in UTM coordinates, Landsat scenes will be clipped 
aoi.bbox <- as.numeric(unlist(strsplit(rciop.getparam("extent"), ",")))
aoi.extent <- extent(aoi.bbox[1], aoi.bbox[3], aoi.bbox[2], aoi.bbox[4])

dest.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

# read the inputs coming from stdin
f <- file("stdin")
open(f)

setwd(TMPDIR)

# create the workspace on GeoServer
workspace.name <- "Test.GeoServer"
CreateGeoServerWorkspace(geoserver, workspace.name)

while(length(ls8.ref <- readLines(f, n=1)) > 0) {

	rciop.log("INFO", paste("processing product", ls8.ref))

	ls8.url <- rciop.casmeta(field="dclite4g:DataSet", url=ls8.ref)$output
	ls8.downloadUrl <- rciop.casmeta(field="dclite4g:onlineResource", url=ls8.ref)$output
	ls8.identifier <- strsplit(rciop.casmeta(field="dc:identifier", url=ls8.ref)$output, ":")[[1]][2]

	# download Landsat 8 product
	rciop.log("INFO", paste("downloading", ls8.url, "to", ls8.identifier))
	rciop.copy (url=ls8.downloadUrl, target=TMPDIR, uncompress=FALSE)

	# extract the compressed Landsat 8 product
	rciop.log("INFO", paste("extracting product", ls8.identifier))
	untar(paste(TMPDIR, "/", ls8.identifier ,".tar.gz", sep=""), exdir = ls8.identifier)

	# read the data
	rciop.log("INFO", paste("Loading", ls8.identifier, "dataset"))
	ls8 <- ReadLandsat8(ls8.identifier, aoi.extent)

	r <- raster()
	coverage.name <- ""
	if (GetOrbitDirection(ls8) == 'A') {
                # ascending direction, get AtSatelliteBrightnessTemperature from TIRS1 
                r <- ToAtSatelliteBrightnessTemperature(ls8, band="tirs1")
                coverage.name <- "Thermal"
    } else {
                # descending direction, get RGB from "swir2", "nir", "green" bands
                r <- ToRGB(ls8, "swir2", "nir", "green")
                coverage.name <- "Colours"
    }
	
	coverage.store <- paste(coverage.name, ls8.identifier, sep="_")               
	rciop.log("INFO", paste0("coverage.store = ",coverage.store))   
    r <- projectRaster(r, crs=dest.proj)
    
    # save the raster to hdsf
    tmp.file <- paste(coverage.store, "tif", sep =".")
    rciop.log("INFO", paste("Save raster to", tmp.file, "file"))
  	writeRaster(r, filename=tmp.file, format="GTiff", overwrite=TRUE)
    res <- rciop.publish(tmp.file, recursive=FALSE, metalink=TRUE)
  	if (res$exit.code==0) { published <- res$output }
	rciop.log("INFO", paste("File", tmp.file, "saved"))

	CreateGeoServerCoverageStore(geoserver,
			                    workspace.name,
            			        coverage.store,
                                TRUE,
                                "GeoTIFF",
                                "file:data/raster.tif")
    
    rciop.log("INFO", "Executing POSTRaster")
    POSTraster(geoserver, workspace.name, coverage.store, r)


	# clean up
	rciop.log("INFO", "Cleaning-up")
	file.remove(tmp.file)
	unlink(paste(TMPDIR, ls8.identifier, sep="/"), recursive=TRUE)

}
