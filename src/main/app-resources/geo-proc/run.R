#!/usr/bin/Rscript --vanilla --slave --quiet

# load rciop library to access the developer cloud sandbox functions
library("rciop")
library("rgeos")
library("stringr")

# load the application package when mvn installed it
library(rGeoServer, lib.loc="/application/share/R/library")
library(rLandsat8, lib.loc="/application/share/R/library")


# get the GeoServer REST access point
geoserver <- rciop.getparam("geoserver")

# get the extent of the area of interest in UTM coordinates, Landsat scenes will be clipped 
aoi.bbox <- as.numeric(unlist(strsplit(rciop.getparam("extent"), ",")))
aoi.extent <- extent(aoi.bbox[1], aoi.bbox[3], aoi.bbox[2], aoi.bbox[4])

country.code <- "ITA"
dest.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

# read the inputs coming from stdin
f <- file("stdin")
open(f)

setwd(TMPDIR)

while(length(ls8.ref <- readLines(f, n=1)) > 0) {
  	# create the country workspace on GeoServer
	CreateGeoServerWorkspace(geoserver, country.code)

	rciop.log("INFO", paste("processing product", ls8.ref))

	ls8.url <- rciop.casmeta(field="dclite4g:DataSet", url=ls8.ref)$output
	ls8.downloadUrl <- rciop.casmeta(field="dclite4g:onlineResource", url=ls8.ref)$output
	ls8.identifier <- strsplit(rciop.casmeta(field="dc:identifier", url=ls8.ref)$output, ":")[[1]][2]

	# download Landsat 8 product
	rciop.log("INFO", paste0("downloading", ls8.url, "to", ls8.identifier))
	rciop.copy (url=ls8.downloadUrl, target=TMPDIR, uncompress=FALSE)

	# extract the compressed Landsat 8 product
	rciop.log("INFO", paste("extracting product", ls8.identifier))
	untar(paste(TMPDIR, "/", ls8.identifier ,".tar.gz", sep=""), exdir = ls8.identifier)

	# read the data
	rciop.log("INFO", paste0("Loading", ls8.identifier, "dataset"))
	ls8 <- ReadLandsat8(ls8.identifier, aoi.extent)

	if (GetOrbitDirection(ls8) == 'A') {
		rciop.log("INFO", "Ascending orbit, saving TIRS1 band")

		# ascending direction, get AtSatelliteBrightnessTemperature from TIRS1 
		bt <- ToAtSatelliteBrightnessTemperature(ls8, band="tirs1")
		r <- projectRaster(bt, crs=dest.proj)

	} else {

		rciop.log("INFO", "Descending orbit, saving RGB image")

		# descending direction, get RGB from "swir2", "nir", "green" bands
		raster.image <- ToRGB(ls8, "swir2", "nir", "green") 
		r <- projectRaster(raster.image, crs=dest.proj)    
	}


	# define the coverage store name
	# coverage.store <- paste("sla", format(as.Date(coverages$start[i]), format="%Y-%m"), sep="_")
	coverage.store <- paste("Etna", ls8$metadata$date_acquired, format="%Y-%m"), sep="_")
	CreateGeoServerCoverageStore(geoserver, 
		                        country.code,
		                        coverage.store,
		                        TRUE,
		                        "GeoTIFF",
		                        "file:data/raster.tif")

	POSTraster(geoserver, country.code, coverage.store, r)


	# publish it
	res <- rciop.publish(ls8.png, recursive=FALSE, metalink=TRUE)
	if (res$exit.code==0) { published <- res$output }

	# publish it
	res <- rciop.publish(ls8.tif, recursive=FALSE, metalink=TRUE)
	if (res$exit.code==0) { published <- res$output }

	# clean up
	rciop.log("INFO", "Cleaning-up")
	unlink(paste(TMPDIR, ls8.identifier, sep="/"), recursive=TRUE)

}
