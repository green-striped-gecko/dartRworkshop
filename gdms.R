library(gdm)

load(system.file("./data/gdm.RData", package="gdm"))# columns 3-7 are soils variables, remainder are climategdmExpData[1:3,]


gdmExpData[1:10,1:5]

names(gdmExpData)


sppTab <- gdmExpData[,c("species", "site", "Lat", "Long")]

envTab <- gdmExpData[,c(2:ncol(gdmExpData))]

gdmTab <-formatsitepair(sppTab, bioFormat=2, XColumn="Long", YColumn="Lat",sppColumn="species", siteColumn="site", predData=envTab)


gdmTab[1:3,]



# environmental raster data
rastFile <-system.file("./extdata/stackedVars.grd", package="gdm")
envRast <-stack(rastFile)
gdmTab.rast <-formatsitepair(sppTab, bioFormat=2, XColumn="Long", YColumn="Lat",sppColumn="species", siteColumn="site", predData=envRast)



gdm.1 <-gdm(gdmTab, geo=T)

library(dartR)
library(PopGenReport)
library(raster)


## Not run: 
data(possums.gl)
library(raster)  #needed for that example
landscape.sim <- readRDS(system.file("extdata","landscape.sim.rdata", package="dartR"))
glc <- gl.genleastcost(x=possums.gl,fric.raster=landscape.sim , 
                       gen.distance = "D", NN=8, pathtype = "leastcost",plotpath = TRUE)


ff <- landscape.sim
values(ff)<- values(ff)*runif(62500)
ff2 <- landscape.sim
values(ff2)<- values(ff2)*runif(62500)


frs <- stack(landscape.sim, ff, ff2)




glc <- gl.genleastcost(x=possums.gl,fric.raster=fric.raster , gen.distance = "D", NN=8, pathtype = "leastcost",plotpath = TRUE)

gl.gdm <- function(x, glc, geo=TRUE, splines=NULL, knots=NULL)
{
  gdis <- data.frame(site=colnames(glc$gen.mat), glc$gen.mat)
  
  #recover coordinates
  fac <- pop(x)
  xys <- apply(x@other$xy, 2, function(q) tapply(q, fac, mean))
  predtab <- data.frame(site=rownames(xys), x=xys[,1], y=xys[,2] )
  costmats <- lapply(glc$cost.mats, function (x) data.frame(site=colnames(xys), x )
  )
  gtab <- formatsitepair(bioData = gdis, bioFormat = 3,XColumn = "x", YColumn = "y", pred=predtab, siteColumn = "site", distPreds = costmats)
  gdmmodell <-gdm(gtab, geo=geo, splines=splines, knots=knots) 
  return(gdmmodell)
}
  

gdis <- data.frame(site=colnames(glc$gen.mat), glc$gen.mat)

predtab <- data.frame(site=names(xs), x=xs, y=ys )

distP1 <- data.frame(site=names(xs), glc$cost.mats[[1]])


rr <- glc$eucl.mat*runif(100)
distP2 <- data.frame(site=names(xs),  rr)
distP3 <- data.frame(site=names(xs),  glc$eucl.mat*runif(100))


gtab <- formatsitepair(bioData = gdis, bioFormat = 3,XColumn = "x", YColumn = "y", pred=predtab, siteColumn = "site", distPreds = list(layer=distP, eucl=distP2, distP3 ))




gdm.1 <-gdm(gtab, geo=F)
summary(gdm.1)
plot(gdm.1)


modTest <- gdm.varImp(gtab, geo=F, nPerm = 50, parallel = T, cores=5)

barplot(sort(modTest[[2]][,1], decreasing=T))

