##################################
### Run Circuitscape (Windows) ###
##################################

# load raster library
library(raster)


gl.runCS <- function(landscape, locs,outpath=tempdir(), plot=FALSE, CS_exe = 'C:/"Program Files"/Circuitscape/cs_run.exe' )
{
# Cost surface
cost <- landscape
# Locs
sites <- locs

# Plot it
if (plot)
{
plot(cost)
points(sites,pch=19,col=2)
}
# Rasterize points using the cost extent

sites <- rasterize(x = sites,y = cost)
# Write rasters to your working directory

writeRaster(sites,file.path(outpath, "sites_rast.asc"),overwrite=TRUE)
writeRaster(cost,file.path(outpath,"cost_rast.asc"),overwrite=TRUE)

# Make an .ini file
CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            "write_cur_maps = 1",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
                  paste(outpath,c("sites_rast.asc",
                                  "cost_rast.asc",
                                  "CS.out"),
                        sep="/")))

# Write it to your working directory
writeLines(CS_ini,file.path(outpath,"myini.ini"))

# Make the CS run cmd
CS_run <- paste(CS_exe, file.path(outpath,"myini.ini"))

# Run the command
system(CS_run)

# Import the effective resistance
CSdis <- as.dist(read.csv(file.path(outpath,"CS_resistances.out"),sep=" ",row.names=1,header=1))

CS_currentmap <- raster(file.path(outpath, "CS_cum_curmap.asc"))

return(list(CSdis=CSdis, map=CS_currentmap))
}





