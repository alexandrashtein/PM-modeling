# This code converts the MAIAC aqua hdf files first into tiff files and finaly into csv file
rm(list=ls())
library(raster)
library(rgdal)
library(rgeos)
library(gdalUtils)
library(plyr)
library(leaflet)
library(magrittr)
library(sp)
library(gdalUtils)

# World
# world = getMap(resolution = "high")
# pol = world[!is.na(world$NAME) & world$NAME == "Switzerland", ]
# pol = spTransform(pol, "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# pol = gBuffer(pol, width = 10000)
# pol = spTransform(pol, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

hdf_to_tif = function(path, n) {
  
  sds = get_subdatasets(path)
  x = strsplit(sds[n], ":")[[1]]
  x = paste0(x[1], ":", x[2], ':"', x[3], '":', x[4], ":", x[5], ":", x[6])
  system(
    paste0(
      "gdal_translate -of GTiff ",
      "\"", 
      x,
      "\" ",
      gsub(".hdf", paste0("_", n, ".tif"), path, fixed = TRUE)
    )
  )
  raster(gsub(".hdf", paste0("_", n, ".tif"), path, fixed = TRUE)) # Convert to raster
  
}

pol=readOGR(dsn="/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/General/Project_border/Project_aoi","Project_border_latlon")
pol = spTransform(pol, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Input Directories
aod_dir = "/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/raw/MAIAC_data_082016"
ref_grid = "/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/raw/MAIAC_data_082016/MAIACLatlon.h03v03.hdf"

# directory aqua old data: N:\Data\Israel\MAIAC_new_10  _2014\h06v07

###################################################################
# STEP 1 - Reading HDF files

# t = "h01v02"
  
# Read 'lon' and 'lat' rasters from static grid
# grid = get_subdatasets(ref_grid)
# lon = grid[2] %>% readGDAL %>% raster
# lat = grid[1] %>% readGDAL %>% raster

sds = get_subdatasets(ref_grid) # Read current file

lon = hdf_to_tif(ref_grid, 2)
lat = hdf_to_tif(ref_grid, 1)

# Creare 'row' and 'col' rasters
row = lon
row[] = rowFromCell(lon, 1:ncell(lon))
col = lon
col[] = colFromCell(lon, 1:ncell(lon))

# Combine to multi-band raster
grid = stack(row, col, lon, lat)
names(grid) = c("row", "col", "lon", "lat")

# Convert to data.frame
grid = as.data.frame(grid)

# Spatial subset
coordinates(grid) = ~ lon + lat
proj4string(grid) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
grid = grid[pol, ]

#####################
# Write grid table

# leaflet() %>% addTiles %>% addPolygons(data = gConvexHull(grid[sample(1:nrow(grid), 1000), ]))
# writeOGR(grid, ".", "grid_h01v02", driver = "ESRI Shapefile")

#####################

# Convert to data.frame
grid = as.data.frame(grid)

for(year in 2015) {

  # Read HDF files list from AOD directory
  setwd(file.path(aod_dir, year))
  
  files = list.files(pattern = "MAIACAAOT.h03v03.*\\.hdf$", recursive = TRUE) # note that MAIACTAOT is for TERRA data and MAIACAAOT is for AQUA data
  
  result = list()
  
  for(f in files) {
  
    # Read data
    sds = get_subdatasets(f)
    # 
    Optical_Depth_047 = hdf_to_tif(f, grep("grid1km:Optical_Depth_047", sds))
    Optical_Depth_055 = hdf_to_tif(f, grep("grid1km:Optical_Depth_055", sds))
    AOT_Uncertainty = hdf_to_tif(f, grep("grid1km:AOT_Uncertainty", sds))
    AOT_QA = hdf_to_tif(f, grep("grid1km:AOT_QA", sds))
    RelAZ=hdf_to_tif(f, grep("grid5km:RelAZ", sds))
    RelAZ=disaggregate(RelAZ, fact = 5)
    
    # Optical_Depth_Land = sds[grepl("Optical_Depth_Land", sds)] %>% readGDAL %>% raster
    # AOT_Uncertainty = sds[grepl("AOT_Uncertainty", sds)] %>% readGDAL %>% raster
    # AOT_QA = sds[grepl("AOT_QA", sds)] %>% readGDAL %>% raster
    row = Optical_Depth_047
    row[] = rowFromCell(Optical_Depth_047, 1:ncell(Optical_Depth_047))
    col = Optical_Depth_047
    col[] = colFromCell(Optical_Depth_047, 1:ncell(Optical_Depth_047))
    r = stack(row, col, Optical_Depth_047, Optical_Depth_055, AOT_Uncertainty, AOT_QA, RelAZ)
    names(r) = c("row", "col", "Optical_Depth_047","Optical_Depth_055", "AOT_Uncertainty", "AOT_QA","RelAZ")
    r = as.data.frame(r)
        
    # Join with 'grid'
    r = join(r, grid, c("row", "col"))
    r = r[!is.na(r$lon) & !is.na(r$lat), ]
    # r = r[indices, ]
    
    # Add filename
    r$date = 
      f %>% 
      strsplit("\\.") %>% 
      sapply("[", 3) %>% 
      substr(1, 7) %>% 
      as.Date(format = "%Y%j")
    
    # Combine results
    result[[f]] = r
  
  }
  
  result = do.call(rbind.fill, result)
  
  write.csv(
    result,
    paste0("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/MAIAC_data_082016/MAIACAAOT_Israel_", year, ".csv"),
    row.names = FALSE
    )

}


###################################################################
# STEP 2 - Processing

# library(plyr)
# 
# aod_dir = "/home/michael/MAIAC_AOD_P024/2011"
# setwd(aod_dir)
# 
# file = 
#   list.files(pattern = paste0("MAIACTAOT_", t, "_.*\\.csv$"))
# 
# dat = read.csv(i, stringsAsFactors = FALSE)
# 
# 
# # dat$row = NULL
# # dat$col = NULL
# 
# cols = c(
#   "MESA_Study_Site", 
#   "row", 
#   "col", 
#   "lon", 
#   "lat", 
#   "date",
#   "Optical_Depth_047", 
#   "AOT_Uncertainty", 
#   "AOT_QA", 
#   "file"
#   )
# dat = dat[, cols]
# 
# y = 
#   i %>% 
#   strsplit("_") %>% 
#   sapply("[", 3) %>% 
#   strsplit("\\.") %>% 
#   sapply("[", 1)
# 
# write.csv(
#   dat, 
#   paste0("MAIACTAOT_Switzerland_final.csv"), 
#   row.names = FALSE
#   )
#     
###################################################################













