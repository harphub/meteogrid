#--------------------------------------#
# Part of R-package geogrid            #
# © Alex Deckmyn                       #
# Released under GPL-3 license         #
#--------------------------------------#

#####################################
### VECTOR FIELDS                 ###
#####################################

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

limage.geofield <- function(field,smooth=FALSE,drawmap=TRUE,
                   maplwd=.5,mapcol='black',
                   map.database="worldmap",...){
  gdomain <- attr(field,"domain")
  glimits <- DomainExtent(gdomain)
  x <- seq(glimits$x0,glimits$x1,length=glimits$nx)
  y <- seq(glimits$y0,glimits$y1,length=glimits$ny)
  z <- field[1:gdomain$nx,1:gdomain$ny]

  limage.default(x=x,y=y,z=z,smooth=smooth,...)

### Iw we add dx/2 at the borders, ini a LatLon this means we shift the meridian?
### in fact not!
### BUT: for a global, VERTEX CENTERED grid, we should add the last

  assign(".Last.domain",gdomain,envir=globalenv())
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
      if(as.numeric(R.Version()$major) < 3)
        .Internal(filledcontour(x,y,z,as.double(levels),col = col ))
      else
        .filled.contour(x,y,z,as.double(levels),col)
    }
    else {
### useRaster is available from v 2.13.0 on. It improves the quality enormously.
### but you may not want to use it if you are exporting to some other output device...
      if(as.numeric(R.Version()$major) + .01*as.numeric(R.Version()$minor) >= 2.13){
        image(x,y,z,xlab="",ylab="",axes=FALSE,xlim=xlim,col=col,
          breaks=levels,asp=asp,useRaster=useRaster,...)
      }
      else{
       image(x,y,z,xlab="",ylab="",axes=FALSE,xlim=xlim,col=col,
          breaks=levels,asp=asp,...)
      }
    }

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

"contour.geofield" <-
function(field,mapcol="black",add=FALSE,drawmap=!add,maplwd=.5,
         map.database="worldmap",...){
  gdomain <- attr(field,"domain")
  glimits <- DomainExtent(gdomain)
  x <- seq(glimits$x0,glimits$x1,length=gdomain$nx)
  y <- seq(glimits$y0,glimits$y1,length=gdomain$ny)

  if(drawmap)
      plot(gdomain,maplwd=maplwd,mapcol=mapcol,add=add,
           drawmap=drawmap,map.database=map.database)
### a future using lattice:
#  contourplot(x, y, field[1:gdomain$nx, 1:gdomain$ny],
#          xlab = "", ylab = "", axes = FALSE, add = ifelse(drawmap,TRUE,add), ...)

contour(x, y, field[1:gdomain$nx, 1:gdomain$ny],
          xlab = "", ylab = "", axes = FALSE, add = ifelse(drawmap,TRUE,add), ...)
}

############################
### SHORTCUTS            ###
############################

iview <- function(x,nlevels=15,color.palette=irainbow,
          title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
          legend=FALSE,breaks=seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),
          length=nlevels),mask=NULL,...){
  if(!is.null(mask)){
    if(is.character(mask)) mask <- eval(parse(text=mask))
    else if (is.expression(mask)) mask <- eval(mask)
    x[eval(expression(mask))] <- NA
  }
  limage(as.geofield(x),color.palette=color.palette,
      plot.title=title(main=title),
      smooth=FALSE,legend=legend,nlevels=nlevels,...)
}

### a new version using Lattice (experimental)
iview2 <- function(x,nlevels=15,color.palette=irainbow,
          title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
          legend=list(space="right",width=1,labels=list(cex=0.8)),
          breaks=seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=nlevels),
          asp=1,mask=NULL,
          drawmap=TRUE,maplwd=.5,mapcol='black',
          map.database="worldmap", ...){
  require(lattice)
  if(!is.null(mask)){
    if(is.character(mask)) mask <- eval(parse(text=mask))
    else if (is.expression(mask)) mask <- eval(mask)
    x[eval(expression(mask))] <- NA
  }
  gdomain <- attributes(x)$domain
  glimits <- DomainExtent(gdomain)
  xc <- seq(glimits$x0,glimits$x1,length=glimits$nx)
  yc <- seq(glimits$y0,glimits$y1,length=glimits$ny)
  grid <- expand.grid(x=xc,y=yc)
  grid$z <- as.vector(x[1:gdomain$nx,1:gdomain$ny])

  if(is.null(breaks)) breaks=pretty(grid$z,nlevels)
  assign(".Last.domain",gdomain,envir=globalenv())


  ppp <- levelplot(z~x*y,data=grid,col.regions=color.palette(2*nlevels),
         main=list(label=title),
         colorkey=legend,cuts=nlevels,
         scales=list(draw=FALSE),xlab="",ylab="",at=breaks,asp=asp,,...)
### problem: the levelplot is only "printed" if it is the last statement of the function...
### and now we should add the map...
  print(ppp)
    if(drawmap) plot(gdomain,add=TRUE,drawmap=TRUE,
         add.dx=TRUE,box=TRUE,maplwd=maplwd,mapcol=mapcol,
         map.database=map.database)

}

cview <- function(x,nlevels=15,
           title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
           mask=NULL,...){
  if(!is.null(mask)){
    if(is.character(mask)) mask <- eval(parse(text=mask))
    else if (is.expression(mask)) mask <- eval(mask)
    x[eval(expression(mask))] <- NA
  }

  contour(as.geofield(x),nlevels=nlevels,main=title,...)
}

fcview <- function(x,nlevels=15,color.palette=irainbow,
           title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
           legend=TRUE,breaks=seq(min(x,na.rm=TRUE),
           max(x,na.rm=TRUE),length=nlevels),mask=NULL,...){
  if(!is.null(mask)){
    if(is.character(mask)) mask=eval(parse(text=mask))
    else if (is.expression(mask)) mask=eval(mask)
    x[eval(expression(mask))]=NA
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

plot.geodomain <- function(x=.Last.domain,add=TRUE,
             maplwd=1,mapcol='black',
             add.dx=TRUE,drawmap=!add, box=drawmap,
             map.database="worldmap",...){
### consistency

  if(add) domain <- .Last.domain
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
    xlim <- c(glimits$x0-glimits$dx/2,glimits$x1+glimits$dx/2)
    ylim <- c(glimits$y0-glimits$dy/2,glimits$y1+glimits$dy/2)
  }

  if(!add){
     x <- seq(glimits$x0,glimits$x1,length=domain$nx)
     y <- seq(glimits$y0,glimits$y1,length=domain$ny)
     image(x,y,array(0,dim=c(domain$nx,domain$ny)),xlab="",ylab="",axes=FALSE,col=0)
     assign(".Last.domain",domain,envir=globalenv())
  }

  if(drawmap){
#    require(maps)
    if(!exists(map.database)) data(list=map.database)
    map <- eval(parse(text=map.database))
### BUG: this simple restriction may leave artificial lines!
#    boundaries <- map[map$x >= glimits$lonlim[1] &
#                    map$x <= glimits$lonlim[2] &
#                    map$y >= glimits$latlim[1] &
#                    map$y <= glimits$latlim[2] ,]
# SO: first call to map.restrict is to reduce the number points to project.
# there may be more subtle/efficient ways to do this, but at least it works.
    boundaries <- map.restrict(map,xlim=glimits$lonlim,ylim=glimits$latlim)
    geo <- project(boundaries, proj = domain$projection,inv = FALSE)
    xyper <- periodicity(domain)
### make sure that no points fall outside the map domain.
#    cat("xlim=",xlim,"\nylim=",ylim,"\n")
    geo=map.restrict(geo,xlim=xlim,ylim=ylim,xper=xyper$xper,yper=xyper$yper)
#    trellis.focus("panel",1,1)
#    panel.lines(geo, col = mapcol, lwd = maplwd,...)
### or use panel.lines ?
    lines(geo, col = mapcol, lwd = maplwd,...)
  }

  if(box){
#    panel.lines(c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),
#          c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),col="black",...)
    lines(c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),
          c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),...)
  }
#  trellis.unfocus()
}

plot.geofield <- function(field,...){
  plot(attr(field,"domain"),...)
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

  if(is.geofield(domain)) domain=attributes(domain)$domain

  glimits <- DomainExtent(domain)
  if (!add.dx) {
    xlim <- c(glimits$x0, glimits$x1)
    ylim <- c(glimits$y0, glimits$y1)
  }
  else {
    xlim <- c(glimits$x0 - glimits$dx/2, glimits$x1 + glimits$dx/2)
    ylim <- c(glimits$y0 - glimits$dy/2, glimits$y1 + glimits$dy/2)
  }


  x <- seq(xlim[1], xlim[2], length = 200)
  y <- seq(ylim[1], ylim[2], length = 200)

  domainframe <- cbind(c(rep(xlim[1],length(y)),x,rep(xlim[2],length(y)),rev(x)),
                         c(y,rep(ylim[2],length(x)),rev(y),rep(ylim[1],length(x)) ) )
  domainlalo <- project(domainframe,proj=domain$projection,inv=TRUE)

  lines(project(domainlalo,proj=.Last.domain$projection),...)

}

############################
### adding point values ####
############################
### this will plot point values on a map
### with a colour function
obsplot <- function(x,y,z,breaks=5,pretty=TRUE,legend.pos=NULL,add=TRUE,domain=.Last.domain,col=irainbow,...){
  if (length(breaks) == 1L & pretty) breaks <- pretty(z,breaks)
  bins <- cut(z,breaks,include.lowest=TRUE,right=FALSE)

  nlev <- length(levels(bins))
  cols <- col(nlev)
  xyp <- project(x,y)
  if(!add) plot.domain(domain,add=FALSE)
  points(xyp,col=cols[as.integer(bins)],pch=19)
  if(!is.null(legend.pos)) legend(x=legend.pos,legend=rev(levels(bins)),fill=rev(cols))
  return(invisible(list(xyp=xyp,z=z,levels=levels(bins),cols=cols)))

}

