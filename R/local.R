##############################
### Some Standard Domains: ###
##############################
"as.geofield" <-
  function (x=array(NA,dim=c(domain$nx,domain$ny)), domain,time = "", info = list())
{
  if(is.geofield(x)) x
  else {
    if(is.geodomain(x)) {
      domain <- x
      x <- array(NA,dim=c(domain$nx,domain$ny))
    }
    else {
### no domain is given, so we start guessing...
      if (missing(domain)) {
        ddd <- dim(x)[1]
        if (ddd == 128 | ddd == 117)
          domain <- DomainMar
        else if (ddd == 108 | ddd == 97)
          domain <- DomainBel
        else if (ddd == 67 )
          domain <- DomainERA10km
        else if (ddd == 181 | ddd == 192)
          domain <- DomainFRBE
        else if (ddd == 240 | ddd == 229)
          domain <- DomainBel2
        else if (ddd == 300 | ddd == 289 | ddd == 320)
          domain <- DomainFra
        else if (ddd == 291)
          domain <- DomainGLAMEPS
        else if (ddd == 486)
          domain <- DomainGLAMEPSv0
        else if (ddd == 646)
          domain <- DomainGLAMEPSv1
        else if (ddd == 601)
          domain <- DomainINCA
        else return("Unknown domain")
      }
    }
    attr(x,"domain") <- domain
    attr(x,"time") <- time
    attr(x,"info") <- info
    class(x) <- c("geofield",class(x))
    x
  }
}
### some useful projections for us Belgians (official Lambert projections):
lambel2008 <- list(proj="lcc",ellps="GRS80",
               lat_1=49.83333333333333,lat_2=51.166666666666667,
               lat_0=50.797815,lon_0=4.35921583333333333,
               y_0=665262.0,x_0=649328.0
               )
lambel1972 <- list(proj="lcc",a=6378388,f=1/297,
                   lat_1=49.833333333333,lat_2=51.1666666666666,
                   lat_0=90,lon_0=4.367487,
                   x_0=150000.013,y_0=5400088.438)

### Domains commonly used in Belgium:
##you can't use project at startup (Dyn library not loaded yet)
##SW=unlist(project(c(360000,350000),proj=lambel2008,inv=TRUE)),
##NE=unlist(project(c(960000,940000),proj=lambel2008,inv=TRUE))),
DomainINCA <- structure(list(projection=lambel2008,nx=601,ny=591,
                             SW=c( 0.491815928394773 , 47.896068330561 ),
                             NE=c( 9.0040364775994 , 53.1788068443894 )),
                        class="geodomain")

DomainMar <- structure(
                       list(projection = list(proj="lcc", lon_0=353, lat_1=31.5605943600000,
                                              lat_2=31.5605943600000,a=6371229.0,es=0.0),
                            nx = 117, ny = 117,
                            SW = c(339.80916592, 18.0573977600000),
                            NE = c(10.3442195400000,  43.3049020977534)),
                       class = "geodomain")

DomainBel.old <-  structure(
                        list(projection = list(proj="lcc", lon_0=4.55361516000000, lat_1=50.5698986500001,
                                              lat_2=50.5698986500001,a=6371229.0,es=0.0),
                             nx = 97, ny = 97,
                             SW = c(0.108210750000000, 47.4739871500001),
                             NE = c(9.60345170000001, 53.4734923027197)),
                        class = "geodomain")

DomainBel7 <-  structure(
                         list(projection = list(proj="lcc", lon_0=3.90000000000003, lat_1=50.7639245100143,
                              lat_2=50.7639245100143,a=6371229.0,es=0),
                              nx = 229, ny = 229,
                              SW = c(-5.8377382093155, 43.1741738778606),
                              NE = c(17.0809205445003, 57.2521782233093)),
                         class = "geodomain")

DomainBel4 <-  structure(
                         list(projection = list(proj="lcc", lon_0=4.55361516, lat_1=50.56989865,
                              lat_2=50.56989865,a=6371229.0,es=0.0),
                              nx = 181, ny = 181,
                              SW = c(359.777917039402,47.2204916197932),
                              NE = c(10.0338376063779,53.6950685547022)),
                         class = "geodomain")


DomainFra <-  structure(
                        list(projection = list(proj="lcc", lon_0=2.57831001000000, lat_1=46.4688491800000,
                                               lat_2=46.4688491800000,a=6371229.0,es=0.0),
                             nx = 289, ny = 289,
                             SW = c(348.16028064, 33.1436856300000),
                             NE = c(25.0576044700000,56.9571605656515)),
                        class = "geodomain")

DomainFRBE <-  structure(
                        list(projection = list(proj="lcc", lon_0=2.57831001000000, lat_1=46.4688491800000,
                                              lat_2=46.4688491800000,a=6371229.0,es=0.0),
                             nx = 181, ny = 181,
                             SW = c(-6.83896855384285,42.6379529283229),
                             NE = c(18.1316332893156 ,57.5314314771777)),
                        class = "geodomain")

DomainFRBE2 <-  structure(
                        list(projection = list(proj="lcc", lon_0=2.57831001000000, lat_1=46.4688491800000,
                                              lat_2=46.4688491800000,a=6371229.0,es=0.0),
                             nx = 289, ny = 349,
                             SW = c(343.343048103504 , 41.6493913129657),
                             NE = c( 30.782850082984 , 70.428851485199)),
                        class = "geodomain")

DomainGLAMEPS <- structure(
                        list(projection = list(proj = "ob_tran",
                                              o_proj = "latlong", o_lat_p = 40, o_lon_p = 0, lon_0 = 22),
                             nx = 291, ny = 260,
                             SW = c(-20.4660686820290, 17.2535716092582),
                             NE = c(80.2596473212195, 75.4164238926903)),
                        class = "geodomain")

DomainGLAMEPSv0 <- structure(
          list(projection = list(proj = "ob_tran",
                                 o_proj = "latlong", o_lat_p = 45,
                                 o_lon_p = 0, lon_0 = 36),
               nx = 486, ny = 378,
               SW = c(-17.7453535169411, 17.9050524844544),
               NE = c(53.654829245849,77.4627171299656)),
          class = "geodomain")

DomainGLAMEPSv1 <- structure(
          list(projection = list(proj = "ob_tran",
                                 o_proj = "latlong", o_lat_p = 46,
                                 o_lon_p = 0, lon_0 = 30),
               nx = 646, ny = 492,
               SW = c( -23.7829009341615 , 15.7785734266855 ),
               NE = c( 81.2229530812028 , 77.6998129394985 )),
          class = "geodomain")

DomainGLAMEPSv2 <- structure(
          list(projection = list(proj = "ob_tran",
                                 o_proj = "latlong", o_lat_p = 46,
                                 o_lon_p = 0, lon_0 = 30),
               nx = 870, ny = 660,
               SW = c( -23.7829009341615 , 15.7785734266855 ),
               NE = c( 84.0654800534012 , 77.55966361171726)),
          class = "geodomain")

### BUG: at install time, Make.domain etc. can't run yet. So this crashes installation.
#DomainGLAMEPSv2.ald  <- Make.domain("lambert",dxdy=c(8900,8900),nxny=c(853,709),
#                            clonlat = c(-3.514364, 55.229520),
#                            reflat=42.8,reflon=28)
#DomainINDRAv1=subgrid(DomainGLAMEPSv1,x1=335,x2=390,y1=170,y2=230)
#DomainINDRAv2=subgrid(DomainGLAMEPSv2,x1=446,x2=521,y1=226,y2=307)

DomainGLAMEPSv1.ald <- structure(
           list(projection = list(proj = "lcc",
                                  lon_0 = 28,lat_1 = 42.8, lat_2 = 42.8,
                                  a = 6371229, es = 0),
                nx = 629, ny=529,
                SW = c(336.268176218338, 15.1916592437107),
                NE = c(85.3445190672328, 78.2862682959375)),
           class = "geodomain")


`DomainGLAMEPS12fa.bad` <-
structure(list(projection = list(proj = "lcc", lon_0 = 35.0, lat_1 = 45,lat_2 = 45,
    a = 6371229, es = 0),
    nx = 501, ny = 421,
    SW = c(342.682814032803, 16.7232476459216),
    NE = c(57.542440516914, 78.3155853407864)),
    class = "geodomain")



DomainERA10km <- structure(
                       list(projection = list(proj="lcc", lon_0=4.55361516, lat_1=50.56989865,
                                              lat_2=50.56989865, a=6371229.0, es=0.0),
                            nx = 67, ny = 67,
                            SW = c(0.163155782953472, 47.5158190985390),
                            NE = c(9.53269211610237, 53.4367017362904)),
                 class = "geodomain")

##########################


