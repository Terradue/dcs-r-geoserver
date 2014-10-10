#!/usr/bin/Rscript --vanilla --slave --quiet

# load rciop library to access the developer cloud sandbox functions
library("rciop")
library("rgeos")
library(stringr)

# load the application package when mvn installed it
library(rGeoServer, lib.loc="/application/share/R/library")


# read the inputs coming from stdin
f <- file("stdin")
open(f)

setwd(TMPDIR)

while(length(ls8.ref <- readLines(f, n=1)) > 0) {
  # to do: write the application
  rciop.log("INFO", paste("processing product", ls8.ref))
  
}
