library(magrittr)
library(XML)
library(raster)
library(rgdal)

name = "israel"
filename = "l7_israel_2000_to_2002_jun_aug_20cl"
code_dir = "/home/michael/Dropbox/BGU/Itai/landsat_israel"
images_dir = "/media/michael/Elements/" %>% paste0(filename)
selected_months = 1:12

#######################################################################
# Download

setwd(code_dir)

dat = filename %>% paste0(".html") %>% htmlParse

files = dat["//a"] %>% sapply(xmlValue)
files = files[grepl("\\.tar.gz$", files)]
files_names = 
  strsplit(files, "/", fixed = TRUE) %>% 
  sapply("[", 6)

# images_dir = "~/Downloads/l8"

setwd(images_dir)

for(i in 1:length(files)) {
  download.file(
    files[i], 
    file.path(images_dir, files_names[i]), 
    method = "wget", 
    quiet = FALSE
  )
}

#######################################################################
# Unzip

setwd(images_dir)

files = list.files(pattern = "\\.tar\\.gz$")

for(i in files) {
  
  untar(i, exdir = file.path(".", gsub("\\.tar\\.gz$", "", i)))
  
}

#######################################################################
# Calculate NDVI

setwd(images_dir)

dirs = list.files()[file.info(list.files())$isdir]

for(i in dirs) {
  
  red = 
    list.files(i, pattern = "_sr_band3\\.tif$", full.names = TRUE) %>% 
    raster
  
  nir = 
    list.files(i, pattern = "_sr_band4\\.tif$", full.names = TRUE) %>% 
    raster
  
  ndvi = (nir - red) / (nir + red)
  
  cfmask = 
    list.files(i, pattern = "_cfmask.tif$", full.names = TRUE) %>% 
    raster
  
  ndvi[cfmask != 0] = NA # Remove non-'clear' pixels
  
  writeRaster(
    ndvi,
    paste0(i, "_ndvi.tif"),
    format = "GTiff",
    overwrite = TRUE
  )
  
  removeTmpFiles(h=0)
  
}

#######################################################################
# Mosaic

setwd(images_dir)

files = list.files(pattern = "_ndvi.tif$")

days = 
  substr(files, 10, 16) %>% 
  as.Date(format = "%Y%j")

months = 
  days %>% 
  format(format = "%m") %>% 
  as.numeric

files = files[months %in% selected_months]

r = lapply(files, raster)

# Mosaic
r$fun = mean
r$na.rm = TRUE
r$progress = "text"
r$filename = paste0(filename, ".tif")
r$format = "GTiff"
r$overwrite = TRUE

do.call(mosaic, r)


#######################################################################
# Crop

setwd(images_dir)

r = filename %>% paste0("_1.tif") %>% raster
aoi = code_dir %>% readOGR("france_hi")
aoi = spTransform(aoi, proj4string(r))

r = mask(r, aoi, progress = "text")
r = trim(r, progress = "text")

writeRaster(r, filename, format = "GTiff", overwrite = TRUE)










