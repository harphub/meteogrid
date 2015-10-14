#--------------------------------------#
# Part of R-package geogrid            #
# Copyright Alex Deckmyn               #
# Released under GPL-3 license         #
#--------------------------------------#

#####################################
### VECTOR FIELDS                 ###
#####################################

.Last.domain <- local({
    ldom <- NULL
    function(new) if(!missing(new)) ldom <<- new else ldom
})


### plotting a (U,V) vectorfield on a domain !
### don't expect too much...
"vecplot"  <- function(U,...) UseMethod("vecplot")

vecplot.geofield <- function(U,V,add=FALSE,aspcorrect=TRUE,
                     drawmap=TRUE,mapcol="black",maplwd=.5,...){

  gdomain <- attr(U,"domain")
  glimits <- DomainExtent(gdomain)
  x <- seq(glimits$x0,glimits$x1,length=gdomain$nx)
  y <- seq(glimits$y0,glimits$y1,length=gdomain$ny)

### asp=1 assumes the co-ordinates of the map are proportional to those of the vector field
### this is WRONG in lat-lon maps! Actually you should rescale by cos(lat)!

#    if ( attr(U,"domain")$projection[1]=="proj=lalo" & (aspcorrect) ){

  ppp <- gdomain$projection
  if(is.list(ppp)) proj <- ppp$proj
  else proj <- substr(ppp[1],5,nchar(ppp[1]))
  if ( aspcorrect & (proj=="lalo" | proj=="rotlalo") ){
    print("LatLon domain: rescaling U components to correct for local aspect ratio U -> U / cos(lat).")
#      aspect=cos(glimits$clonlat[2])
#      print(paste("Using aspect ratio at domain center:",aspect))
    localaspect <- cos(DomainPoints(U,type="lalo")$lat * pi/180)
### get the speed correct!
    vecnorm <- sqrt(U^2 + V^2)
    U <- U / localaspect
### Correct the norm for the map factor:
    newnorm <- sqrt(U^2 + V^2)
    U <- U*vecnorm/newnorm
    V <- V*vecnorm/newnorm
  }

  if (drawmap)
    plot(gdomain,maplwd=maplwd,mapcol=mapcol,add=add,box=!add,add.dx=TRUE)

  vecplot(U=U[1:gdomain$nx,1:gdomain$ny],
          V=V[1:gdomain$nx,1:gdomain$ny],x=x,y=y,add=TRUE,...)

}

"vecplot.default" <-
function(U,V,x=1:dim(U)[1],y=1:dim(U)[2],thinx=10,thiny=thinx,aspect=1,
         add=FALSE,arrowcolor=1,arrowsize=.03,maxscale=.9,lwd=.5,rescale,xlab="",ylab="",...)
{
### The vectors are plotted as "arrows" (line elements!) in the current co-ordinate system.
### warning("VECPLOT is quite primitive. Use it only for a rough idea of vector fields.")

### The (global) aspect ratio determines the rescaling of the axes. Default is 1.
#  if (missing(aspect)) {
#    aspect <- 1
#    print("Setting aspect=1 (equal units for both axes).")
#  }
  U <- U / aspect

### thinning : not all domain points plotted.
  filterx <- (( (1:length(x) - 1) %% thinx) == 0)
  filtery <- (( (1:length(y) - 1) %% thiny) == 0)
  filterxy <- as.logical(filterx %*% t(filtery))

### the selected points are thus:
  xcoord <- x[as.logical(filterx)]
  ycoord <- y[as.logical(filtery)]

### create a list of the x and y coordinates of all the points in the (filtered) domain
  x <- rep(xcoord,length(ycoord))
  y <- rep(ycoord,rep(length(xcoord),length(ycoord)))
### you don't want arrows to overlap too much
### I simply ensure that u and v are smaller than .9 (maxscale) times the (thinned) domain distance
### can be overruled by rescale
  if(missing(rescale)){
    dx <- xcoord[2]-xcoord[1]
    dy <- ycoord[2]-ycoord[1]
    scale <-  maxscale * min( dx/max(abs(U[filterxy]),na.rm=TRUE) , dy/(max(abs(V[filterxy]),na.rm=TRUE)),na.rm=TRUE )
  }
  else
    scale <- rescale

### The vector components are
### in fact nothing else than the end points of arrows
### v is scaled according to the aspect ratio (asp) of the axes.
### This is quite "synthetic" and not really suitable
### but I see no other way.
  u <- x + scale * as.vector(U[filterxy])
  v <- y + scale * as.vector(V[filterxy])

### create a plot area
  if(!add) plot(c(x[1],x[length(x)]),c(y[1],y[length(y)]),col=0,axes=FALSE,xlab=xlab,ylab=ylab,...)
### at last, the vectors...
### avoid a list of warnings when wind speed is 0...
  suppressWarnings(arrows(x,y,u,v,length=arrowsize,col=arrowcolor,lwd=lwd))
  if(!add) box()
}
###################################
### IMAGE and FILLED.CONTOUR merged
### with legend added
####################################

limage  <- function(x,...) UseMethod("limage")

limage.geofield <- function(x,smooth=FALSE,drawmap=TRUE,
                   maplwd=.5,mapcol='black',
                   map.database="world",...){
  gdomain <- attr(x,"domain")
  glimits <- DomainExtent(gdomain)

  limage.default(x=seq(glimits$x0,glimits$x1,length=glimits$nx),
                 y=seq(glimits$y0,glimits$y1,length=glimits$ny),
                 z=x[1:gdomain$nx,1:gdomain$ny],
                 smooth=smooth,...)

### If we add dx/2 at the borders, in a LatLon this means we shift the meridian?
### in fact not!
### BUT: for a global, VERTEX CENTERED grid, we should add the last

  .Last.domain(gdomain)
  if(drawmap)
    plot(gdomain,add=TRUE,drawmap=TRUE,
         add.dx=!smooth,box=TRUE,maplwd=maplwd,mapcol=mapcol,
         map.database=map.database)
}

limage.default <-
  function(x=1:dim(z)[1],y=1:dim(z)[2],z,smooth=FALSE,
           nlevels=15,levels=pretty(zlim,nlevels),
           color.palette=colorRampPalette(c("blue","white","red")),
           col=color.palette(length(levels)-1),
           legend=FALSE,legend.cex=.7,
           legend.width=1/12,legend.sep=c(1/4,1/2),
           legend.skip=1,legend.digits=5,
           plot.title,title.adjust=TRUE,
           legend.title,
           asp=1, useRaster=TRUE, ...){
    zlim <- range(z,finite=TRUE)
    ncol <- length(levels)-1
    if(smooth){
      xlim <- range(x)
      ylim <- range(y)
    }
    else {
      dx <- diff(x)[1]
      dy <- diff(y)[1]
      xlim <- c(min(x)-dx/2,max(x)+dx/2)
      ylim <- c(min(y)-dy/2,max(y)+dy/2)
    }
    mapxlim <- xlim
    mapylim <- ylim

    if(legend)xlim[2] <- xlim[1]+(xlim[2]-xlim[1])/(1-legend.width)

    nlevels <- length(levels)
    if(smooth){
      plot.new()
      plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp=asp)
      .filled.contour(x,y,z,as.double(levels),col)
    }
    else {
### useRaster is available from v 2.13.0 on. It improves the quality enormously.
### but you may not want to use it if you are exporting to some other output device...
      image(x, y, z, xlab="", ylab="", axes=FALSE, xlim=xlim, col=col,
            breaks=levels, asp=asp, useRaster=useRaster, ...)
    }
### The legend is drawn inside the same plot area, so using the same co-ordinate space as the map.
### That may seem weird, but it is quite effective, fast and easy to combine several plots.
    if(legend){
      legendlevels <- levels
#      if(legendlevels[1]<min(z,na.rm=TRUE))legendlevels[1]=min(z,na.rm=TRUE)
#      if(legendlevels[(ncol+1)]>max(z,na.rm=TRUE))
#            legendlevels[(ncol+1)]=max(z,na.rm=TRUE)
      legw <- legend.width*(mapxlim[2]-mapxlim[1])

      legxlim <- mapxlim[2]+legw*legend.sep
      ybreaks <- seq(ylim[1],ylim[2],length=ncol+1)
      rect(xleft=rep(legxlim[1],ncol),xright=rep(legxlim[2],ncol),
           ybottom=ybreaks[1:ncol],ytop=ybreaks[2:(ncol+1)],
           col=col)

      text.x <- legxlim[2]+(xlim[2]-legxlim[2])/6
      text.y <- c(ylim[1],ybreaks[2:ncol],ylim[2])

      text(x=text.x,y=ylim[1],labels=format(min(legendlevels),digits=5),
           cex=legend.cex,adj=c(0,0))
      y.index <- seq(1+legend.skip,ncol,by=legend.skip)
      if(ncol+1-y.index[length(y.index)]<legend.skip/2)
        y.index=y.index[-length(y.index)]
      text(x=text.x,y=text.y[y.index],
           labels=format(legendlevels[y.index],
           digits=legend.digits),
           cex=legend.cex,adj=c(0,.5))
      text(x=text.x,y=ylim[2],labels=format(max(legendlevels),digits=5),
           cex=legend.cex,adj=c(0,1))
      if(!missing(legend.title)){
        oadj <- par("adj")
        nadj <- oadj+oadj*(1-legend.width)
        par(adj=nadj)
        legend.title
        par(adj = oadj)
      }

    }
    if(!missing(plot.title)){
      if(legend & title.adjust){
        oadj <- par("adj")
        nadj <- oadj*(1-legend.width)
        par(adj=nadj)
      }
      plot.title
      if (legend & title.adjust) par(adj = oadj)
    }
  }


########################
### CONTOUR          ###
########################

cview <- function(x,nlevels=15,
           title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
           mask=NULL,mapcol="black",add=FALSE,drawmap=!add,maplwd=.5,
           map.database="world",...){
  if(!is.null(mask)){
    if(is.character(mask)) mask <- eval(parse(text=mask))
    else if (is.expression(mask)) mask <- eval(mask)
    x[eval(expression(mask))] <- NA
  }

  gdomain <- attr(x,"domain")
  glimits <- DomainExtent(gdomain)

  if (drawmap) plot(gdomain,maplwd=maplwd,mapcol=mapcol,add=add,
                    drawmap=drawmap,map.database=map.database)

  contour(x=seq(glimits$x0,glimits$x1,length=gdomain$nx),
          y=seq(glimits$y0,glimits$y1,length=gdomain$ny),
          z=x[1:gdomain$nx, 1:gdomain$ny],
          xlab = "", ylab = "", axes = FALSE, add = ifelse(drawmap,TRUE,add), ...)
}

############################
### SHORTCUTS            ###
############################

iview <- function(x,nlevels=15,color.palette=irainbow,
          title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
          legend=FALSE,mask=NULL,...){
  if(!is.null(mask)){
    if(is.character(mask)) mask <- eval(parse(text=mask))
    else if (is.expression(mask)) mask <- eval(mask)
    x[eval(expression(mask))] <- NA
  }
  limage(as.geofield(x),color.palette=color.palette,
      plot.title=title(main=title),
      smooth=FALSE,legend=legend,nlevels=nlevels,...)
}

fcview <- function(x,nlevels=15,color.palette=irainbow,
           title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
           legend=TRUE,mask=NULL,...){
  if(!is.null(mask)){
    if(is.character(mask)) mask <- eval(parse(text=mask))
    else if (is.expression(mask)) mask <- eval(mask)
    x[eval(expression(mask))] <- NA
  }
  limage(as.geofield(x),color.palette=color.palette,plot.title=title(main=title),
          smooth=TRUE,legend=legend,nlevels=nlevels,...)
}

vview <- function(U,V,...){
  vecplot(as.geofield(U),as.geofield(V),...)
}

###############################
### PLOTTING A MAP OR FRAME ###
###############################

plot.geodomain <- function(x=.Last.domain(),
             add=!is.null(.Last.domain()),
             maplwd=1,mapcol='black',
             add.dx=TRUE,drawmap=!add, box=drawmap,
             map.database="world",...){

### consistency
  if (add) domain <- .Last.domain()
  else {
    domain <- x
    add.dx <- TRUE
  }

### for backward compatibility
  glimits <- DomainExtent(domain)

  if(!add.dx){
    xlim <- c(glimits$x0,glimits$x1)
    ylim <- c(glimits$y0,glimits$y1)
  }
  else {
    xlim <- c(glimits$x0,glimits$x1) + glimits$dx*c(-1,1)/2
    ylim <- c(glimits$y0,glimits$y1) + glimits$dy*c(-1,1)/2

  }

  if(!add){
     x <- seq(glimits$x0,glimits$x1,length=domain$nx)
     y <- seq(glimits$y0,glimits$y1,length=domain$ny)
     image(x,y,array(0,dim=c(domain$nx,domain$ny)),xlab="",ylab="",axes=FALSE,col=0)
     .Last.domain(domain)
  }

  if(drawmap){
    boundaries <- map(database=map.database,xlim=glimits$lonlim,ylim=glimits$latlim,plot=FALSE)
    geo <- project(boundaries, proj = domain$projection,inv = FALSE)
    xyper <- periodicity(domain)
### make sure that no points fall outside the map domain.
    geo <- map.restrict(geo,xlim=xlim,ylim=ylim,xperiod=xyper$xper,yperiod=xyper$yper)
    lines(geo, col = mapcol, lwd = maplwd,...)
  }

  if(box){
#    panel.lines(c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),
#          c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),col="black",...)
    lines(c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),
          c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),...)
  }
}

plot.geofield <- function(x,...){
  plot(attr(x,"domain"),...)
}
######################
### periodicity depends on the projection
### may be 360, 2*pi or ~40.000 in Mercator !
### a non-periodic LatLon domain still has a "period" of 360 !

map.restrict1 <- function(bx,by,xlim,xperiod=NA_real_,xfrac=0.5){
  lx <- length(bx)

  result <- .C("maprestrict1",bx=as.double(bx),by=as.double(by),lx=as.integer(lx),
            x0=as.double(xlim[1]),x1=as.double(xlim[2]),
            nx=double(2*lx),ny=double(2*lx),
            newlen=integer(1),xperiod=as.numeric(xperiod),
            xfrac=as.numeric(xfrac),NAOK=TRUE)
  data.frame(x=result$nx[2:result$newlen],y=result$ny[2:result$newlen])
}

map.restrict <- function(bxy,xlim,ylim,xperiod=NA_real_,xfrac=0.5,yperiod=NA_real_,yfrac=NA_real_){
  bxy <- map.restrict1(bxy$x,bxy$y,xlim,xperiod,xfrac)
  byx <- map.restrict1(bxy$y, bxy$x, ylim, yperiod, yfrac)
  list(x = byx$y, y = byx$x)
}

#####################################################
### plotting the frame of a domain on another map ###
#####################################################

domainbox <-
  function (domain , add.dx = TRUE, ...)
{
  if (is.null(.Last.domain())) stop("There is no image yet to add the domainbox to.")
  if (is.geofield(domain)) domain <- attributes(domain)$domain

  glimits <- DomainExtent(domain)
  if (!add.dx) {
    xlim <- c(glimits$x0, glimits$x1)
    ylim <- c(glimits$y0, glimits$y1)
  }
  else {
    xlim <- c(glimits$x0, glimits$x1) + glimits$dx*c(-1,1)/2
    ylim <- c(glimits$y0, glimits$y1) + glimits$dy*c(-1,1)/2

  }


  x <- seq(xlim[1], xlim[2], length = 200)
  y <- seq(ylim[1], ylim[2], length = 200)

  domainframe <- cbind(c(rep(xlim[1],length(y)),x,rep(xlim[2],length(y)),rev(x)),
                         c(y,rep(ylim[2],length(x)),rev(y),rep(ylim[1],length(x)) ) )
  domainlalo <- project(domainframe,proj=domain$projection,inv=TRUE)
  lines(project(domainlalo,proj=.Last.domain()$projection),...)

}

############################
### adding point values ####
############################
### this will plot point values on a map
### with a colour function
obsplot <- function(x,y,z,breaks=5,pretty=TRUE,legend.pos=NULL,
                    add=TRUE,domain=.Last.domain(),col=irainbow,...){
  if (missing(y)) {
    if(ncol(x)>1) {
      y <- x[,2]
      x <- x[,1]
    } else {
     y <- x[2]
     x <- x[1]
    }
  }

  if (length(breaks) == 1L & pretty) breaks <- pretty(z,breaks)
  bins <- cut(z,breaks,include.lowest=TRUE,right=FALSE)

  nlev <- length(levels(bins))
  cols <- col(nlev)
  xyp <- project(x,y)
  if (!add) plot.geodomain(domain,add=FALSE)
  points(xyp,col=cols[as.integer(bins)],pch=19)
  if (!is.null(legend.pos)) legend(x=legend.pos,legend=rev(levels(bins)),fill=rev(cols))
  return(invisible(list(xyp=xyp,z=z,levels=levels(bins),cols=cols)))

}

############################################
### Add latitude/longitude grid to a map ###
############################################

DrawLatLon <- function(nx=9, ny=9, lines=TRUE, labels=TRUE, cex=1, col="grey",
                       lty=2, font=2, lab.x=2, lab.y=2, ...) {
  if (is.null(.Last.domain())) stop("Sorry, no projection has been defined.")
  glimits <- DomainExtent(.Last.domain())
  xlim <- glimits$lonlim
  ylim <- glimits$latlim

  x <- pretty(xlim, nx)
  lonlist <- x[which(x >= xlim[1] & x <= xlim[2])]

  y <- pretty(ylim, ny)
  latlist <- y[which(y >= ylim[1] & y <= ylim[2])]

  if (lines) {
    lonlines <- expand.grid(y=c(seq(ylim[1],ylim[2],len=100),NA),x=lonlist)
    latlines <- expand.grid(x=c(seq(xlim[1],xlim[2],len=100),NA),y=latlist)
    lalolines <- rbind(latlines,lonlines)
    plines <- project(lalolines,proj=.Last.domain()$projection)
    plines <- map.restrict(plines, c(glimits$x0,glimits$x1),
                                   c(glimits$y0,glimits$y1),
                                   xperiod=periodicity(.Last.domain())$xper)
    lines(plines, col=col, lty=lty)
  }
  if (labels) {
# labels on axis
    NN <- 100
    DX <- diff(glimits$lonlim)/(NN-1) # not a good criterion...
    DY <- diff(glimits$latlim)/(NN-1)

    tx <- data.frame(x=seq(glimits$x0,glimits$x1,length=NN), y=rep(glimits$y0, NN))
    ptx <- project(tx,proj=.Last.domain()$projection,inv=TRUE)
    zx <- vapply(lonlist,function(ll) which.min(abs(ptx$x-ll)),1)
    at.x <- ifelse(abs(ptx$x[zx]-lonlist) < DX, tx$x[zx], NA)

    ty <- data.frame(x=rep(glimits$x0, NN), y=seq(glimits$y0,glimits$y1,length=NN))
    pty <- project(ty,proj=.Last.domain()$projection,inv=TRUE)
    zy <- vapply(latlist,function(ll) which.min(abs(pty$y-ll)),1)
    at.y <- ifelse(abs(pty$y[zy]-latlist) < DY, ty$y[zy], NA)

    axis(1,at=at.x, labels=format(lonlist), tick=!lines, line=-0.5, col.axis=col)
    axis(2,at=at.y, labels=format(latlist), tick=!lines, line=-0.5, col.axis=col)

  }
}
