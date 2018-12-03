
create.pops <- function(n, mindist, landscape, plot=TRUE)
{  
  minx <- raster::extent(landscape)@xmin #get the min and max  coordinates
  miny <- raster::extent(landscape)@ymin #coordinates of the landscape
  maxx <- raster::extent(landscape)@xmax
  maxy <- raster::extent(landscape)@ymax
  
  cc <- 1
  coords <- data.frame(lx=NA, ly=NA)
  while (cc<= n )  #repeat until you have found n locations
  {
    draw=FALSE
    while (draw==FALSE)
    {
      x <- runif(1,minx,maxx)
      y <- runif(1,miny,maxy)
      if (landscape[cellFromXY(landscape,c(x,y) )]==1)  draw=TRUE 
      #check if in the habitat
    }
    
    coords[cc,] <- c(x,y)
    
    if (nrow(coords)>1) d <- min(dist(coords)) else d <- mindist+1 
    
    if (d > mindist) cc <- cc+1  
    #take the location only if distance is larger than mindist
  }
  if (plot==TRUE) 
  {
    plot(landscape)  
    points(coords, pch=16)
  }
  return( as.matrix( coords))
}

pops2gl <- function(simpops, locs=NULL, n.cov=3)
{
simpops <- simpops
locs <- locs
n.cov=n.cov
gi <- pops2genind(simpops, locs, n.cov)
gl <- gi2gl(gi)
}



create.resistance <- function(nx=50, ny=50, p=0.5, A=0.5, resVal=10, plotres=TRUE) {
nx=nx
ny=ny
tempmask<-secr::make.mask(nx=nx,ny=ny,spacing=1)
r <- secr::raster(randomHabitat(tempmask, p = p, A = A))
#set non-habitat to friction values of 10
values(r)[is.na(values(r))==T]<- resVal
if (plotres) plot(r)
return(r)
}

ggmap_to_raster <- function(x){
  if(!inherits(x, "ggmap")) ps_error("x must be a ggmap object (e.g. from ps_bbox_ggmap)")
  
  map_bbox <- attr(x, 'bb')
  .extent <- raster::extent(as.numeric(map_bbox[c(2,4,1,3)]))
  my_map <- raster::raster(.extent, nrow= nrow(x), ncol = ncol(x))
  rgb_cols <- setNames(as.data.frame(t(col2rgb(x))), c('red','green','blue'))
  red <- my_map
  raster::values(red) <- rgb_cols[['red']]
  green <- my_map
  raster::values(green) <- rgb_cols[['green']]
  blue <- my_map
  raster::values(blue) <- rgb_cols[['blue']]
  ras <- raster::stack(red,green,blue)
  ras
}


