
#elevation data of the world:
library(raster)

ts <- raster("./data/TreeScrub.tif")
res <- raster(".data/srtm_72_20.tif")

respro <- project

# use poss.nz and its projection NZGD2000

xmin <- min(poss.nz@other$xy$x-100)
xmax <- max(poss.nz@other$xy$x+100)

ymin <- min(poss.nz@other$xy$y-100)
ymax <- max(poss.nz@other$xy$y+100)

res.temp <- raster( xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, resolution=100, crs=proj4string(ts))

res.proj <- projectRaster(res, crs = proj4string(ts), res = 1000)

res.proj.clip <- crop(res.proj, res.temp)


