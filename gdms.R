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




glc <- gl.genleastcost(x=possums.gl,fric.raster=frs , gen.distance = "D", NN=8, pathtype = "leastcost",plotpath = TRUE)

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





library(BEDASSLE)

data(HGDP.bedassle.data)
data(mcmc
     
     poss.counts <- apply(as.matrix(possums.gl), 2, function(x) tapply(x, pop(possums.gl), mean, na.rm=T))
     
     poss.sample_sizes <- apply(as.matrix(possums.gl), 2, function(x) tapply(!is.na(x), pop(possums.gl), sum, na.rm=T))
     
     poss.D <- glc$eucl.mat
     poss.E = glc$cost.mats
     poss.k = nPop(possums.gl)
     poss.loci = nLoc(possums.gl)
     
     
     mcmc.operators$delta  <- 3
     
mm <- MCMC(
  counts = poss.counts,
  sample_sizes = poss.sample_sizes,
  D = poss.D,
  E = poss.E,
  k = poss.k,
  loci = poss.loci,
  delta = mcmc.operators$delta,
  aD_stp = mcmc.operators$aD_stp,
  aE_stp = mcmc.operators$aE_stp,
  a2_stp = mcmc.operators$a2_stp,
  thetas_stp = mcmc.operators$thetas_stp,
  mu_stp = mcmc.operators$mu_stp,
  ngen = 1000,
  printfreq = mcmc.operators$printfreq,
  savefreq = mcmc.operators$savefreq,
  samplefreq = mcmc.operators$samplefreq,
  directory = NULL,
  prefix = "poss_",
  continue = FALSE,
  continuing.params = NULL)     
 ho


  
MCMC.outpout <- "poss_MCMC_output1.Robj"

plot_all_acceptance_rates(MCMC.outpout)

plot_all_trace("example_MCMC_output1.Robj", percent.burnin = 0, thinning = 1, population.names = NULL)

plot_all_marginals("example_MCMC_output1.Robj", percent.burnin = 0, thinning = 1,population.names = NULL)



