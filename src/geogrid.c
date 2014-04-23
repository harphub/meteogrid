#include<stdio.h>
#include<math.h>
#include<string.h>
#include "R.h"
void maprestrict1(double* bx,double* by,int* lx,
               double* x0,double* x1,
               double* nx,double* ny,int* newlength,
               double* xperiod, double* xfrac);
/* A routine to make sure that a given map does not extend beyond the borders of the domain in the x direction */
/* because when drawing the lines, you may go outside the domain. */
/*
    INPUT:   bx,by,lx=length
             x0,x1 : boundaries
    OUTPUT : nx,ny,newlength
             The length of the date may change! nx and ny must be long enough to accomodate this.

AUTHOR: Alex Deckmyn

*/
/*
#define DEBUG
*/

void maprestrict1(double *bx,double *by,int* lx,
               double *x0,double *x1,
               double *nx,double *ny,int* newlength,
               double *xperiod, double* xfrac)
{
  int i,j,npos,NPOSMAX;
  double rx,ry ;

#ifdef DEBUG
  Rprintf("NA= %lf\n",NA_REAL);
  if(!ISNA(*xperiod)){
  Rprintf("map.restrict: xperiod= %lf, xfrac=%lf \n",*xperiod,*xfrac);
  }
#endif
  /* if it is not explicitely defined,
     we only consider jumps that are larger than the period
     of the X domain. In other words, after the corrections
     of periodicity, you may be sure that you will never skip.
  */
  if(!ISNA(*xperiod) && ISNA(*xfrac) ){
    Rprintf("setting periodicity jump fraction to 1.\n");
    *xfrac = 1. ;
  }
  if(!ISNA(*xperiod) && *xperiod < *x1 - *x0) {
    Rf_error("INCONSISTENCY: X limits larger than periodicity of domain.\n");
  }
  NPOSMAX=2*(*lx);
  nx[0]=ny[0]=NA_REAL;
  npos=1;
  for(i=0;i<*lx;i++){
    if(!ISNA(bx[i]) && !ISNA(*xperiod)) {
      if(bx[i] < *x0) {
        for(j=1 ; bx[i]+j* *xperiod < *x0 ; j++){ }
        if( bx[i] + j* *xperiod <= *x1) bx[i]=bx[i]+j* *xperiod ;
      }
      if(bx[i] > *x1) {
        for(j=1;bx[i]-j* *xperiod > *x1;j++){}
        if( bx[i] - j* *xperiod >= *x0) bx[i]=bx[i]-j* *xperiod ;
      }

    }
    if(ISNA(bx[i])) { /* it is a NA entry that separates 2 segments */
      nx[npos]=NA_REAL;
      ny[npos++]=NA_REAL;
      if(npos>=NPOSMAX) Rf_error("overflow\n");
    }
/* the following should also be done if the point is outside the domain */
    else if(  !ISNA(nx[npos-1]) && !ISNA(*xperiod) &&
              abs(bx[i]-nx[npos-1])>= *xperiod * *xfrac) {
      /* HERE WE SUPPOSE THE POINT SKIPPED OVER A PERIODIC DOMAIN */
      /* So we have to split the line ! */
#ifdef DEBUG
      Rprintf("Splitting at %d %lf %lf -> %lf %lf \n",npos,nx[npos-1],ny[npos-1],bx[i],by[i]);
#endif
      if(npos >= NPOSMAX-2) Rf_error("overflow\n");
      if (nx[npos-1] < bx[i]) { /* split at x0 */
        ry = by[i] + (by[i]-ny[npos-1])/(*xperiod-bx[i]+nx[npos-1])*(bx[i]-*x1);

        nx[npos]   = *x0;
        ny[npos++] = ry;
        nx[npos]   = NA_REAL;
        ny[npos++] = NA_REAL;
        if(bx[i] <= *x1){
          if(npos >= NPOSMAX-2) Rf_error("overflow\n");
          nx[npos]   = *x1;
          ny[npos++] = ry;
          nx[npos]=bx[i];
          ny[npos++]=by[i];
        }
      }
      else { /* split at x1 */
        ry = by[i] + (by[i]-ny[npos-1])/(*xperiod+bx[i]-nx[npos-1])*(*x0-bx[i]);
        nx[npos]   = *x1;
        ny[npos++] = ry;
        nx[npos]   = NA_REAL;
        ny[npos++] = NA_REAL;
        if(bx[i] >= *x0){
          if(npos >= NPOSMAX-2) Rf_error("overflow\n");
          nx[npos]   = *x0;
          ny[npos++] = ry;
          nx[npos]=bx[i];
          ny[npos++]=by[i];
        }
      }
#ifdef DEBUG
      Rprintf("--------------- md-Y= %lf \n",ry) ;
#endif
    }
    else if (bx[i]<=*x0) { /* the point is outside the domain on the left */
      if ( !ISNA(nx[npos-1]) ) {/* if the previous entry was NA :
                                   first point of a segment, do nothing yet*/
        nx[npos]=*x0;
        ny[npos++]=by[i-1]+(by[i]-by[i-1])/(bx[i]-bx[i-1])*(*x0-bx[i-1]);
        if(npos>=NPOSMAX) Rf_error("overflow\n");
        nx[npos]=NA_REAL;
        ny[npos++]=NA_REAL;
        if(npos>=NPOSMAX) Rf_error("overflow\n");
      }
      if (i<*lx-1) { /* if the next entry is in the domain, find border point */
        if (!ISNA(bx[i+1]) && bx[i+1]>*x0 && bx[i+1]<*x1) {
          if (!ISNA(nx[npos-1])){
            nx[npos]=NA_REAL;
            ny[npos++]=NA_REAL;
            if(npos>=NPOSMAX) Rf_error("overflow\n");
          }
          if( !ISNA(*xperiod) && abs(bx[i]-bx[i+1]) > *xperiod * *xfrac ){
#ifdef DEBUG
            Rprintf("Pre-jump???\n");
#endif
            nx[npos]=*x1;
            ny[npos++] = by[i] + (by[i+1]-by[i])/(bx[i+1]-bx[i]- *xperiod)*(*x1-bx[i]-*xperiod);
          }
          else {
            nx[npos]=*x0;
            ny[npos++]=by[i+1]+(by[i]-by[i+1])/(bx[i]-bx[i+1])*(*x0-bx[i+1]);
          }
          if(npos>=NPOSMAX) Rf_error("overflow\n");
        }
      }
    }
    else if (bx[i]>=*x1) { /* the point is outside the domain on the right */
      if( !ISNA(nx[npos-1])) {
        nx[npos]=*x1;
        ny[npos++]=by[i-1]+(by[i]-by[i-1])/(bx[i]-bx[i-1])*(*x1-bx[i-1]);
        if(npos>=NPOSMAX) Rf_error("overflow\n");
        nx[npos]=NA_REAL;
        ny[npos++]=NA_REAL;
        if(npos>=NPOSMAX) Rf_error("overflow\n");
      }
      if (i<*lx-1) {
        if(!ISNA(bx[i+1]) && bx[i+1]>*x0 && bx[i+1]<*x1 ) {
          if (!ISNA(nx[(npos-1)])){
            nx[npos]=NA_REAL;
            ny[npos++]=NA_REAL;
            if(npos>=NPOSMAX) Rf_error("overflow\n");
          }
          if( !ISNA(*xperiod) && abs(bx[i]-bx[i+1]) > *xperiod * *xfrac ){
#ifdef DEBUG
            Rprintf("Pre-jump???\n");
#endif
            nx[npos]=*x0;
            ny[npos++] = by[i] + (by[i+1]-by[i])/(bx[i+1]-bx[i]+ *xperiod)*(*x0 - bx[i]+ *xperiod);
          }
          else {
            nx[npos]=*x1;
            ny[npos++]=by[i+1]+(by[i]-by[i+1])/(bx[i]-bx[i+1])*(*x1-bx[i+1]);
          }
          if(npos>=NPOSMAX) Rf_error("overflow\n");
          if(npos>=NPOSMAX) Rf_error("overflow\n");
        }
      }
    }
    else { /* a point inside the domain */
      if(npos>=NPOSMAX-1) Rf_error("overflow\n");
      nx[npos]   = bx[i];
      ny[npos++] = by[i];
    }
  }
  *newlength=npos-ISNA(nx[npos-1]);
}
/***********************************/

/* Interpolation for (reduced) gaussian grids */
/* given as a set of lat/lon values and NLON */
/* what we calculate is a set of weights (fractional indices) */

void gridindex(double* inlat, double* inlon, int* nx, double* gridlat, int* nlon,
               int* nlat, double* iout1, double* iout2, double* jout){
  int pp,j;
  double dlon1,dlon2;

/* calculate the longitudes for every latitude */
/* basically just devide 2*pi by Nlon */
  for(pp=0;pp < *nx;pp++){
    j=0;
    while( (j<*nlat) && (gridlat[j]< inlat[pp]) ){j++;};
/* TODO:  what if j=0? */
/* VERY important for Arpege files: the "pole" is Paris.*/
/* polar region: " flag" by giving an index outside the region 1:nlat */
    if (j==0) jout[pp]=0;
    else if( j==*nlat) jout[pp]=*nlat+1 ;
/* notice: In fact 1+(j-1): we want index to start at 1 */
    else jout[pp] = j + (inlat[pp]-gridlat[j-1])/(gridlat[j]-gridlat[j-1]);
    if (j>0) {
      dlon1=360./nlon[j-1];
      iout1[pp] = 1+inlon[pp]/dlon1;
/*    Rprintf("        i1=%lf : %d, d=%lf %lf \n",iout1[pp],nlon[j-1],dlon1,iout1[pp]*dlon1);
*/
    }
    else iout1[pp]=NA_REAL;
    if (j<*nlat) {
      dlon2=360./nlon[j];
      iout2[pp] = 1+inlon[pp]/dlon2;
/*    Rprintf("        i2=%lf : %d,d=%lf %lf\n",iout2[pp],nlon[j],dlon2,iout2[pp]*dlon2);
*/
    }
    else iout2[pp]=NA_REAL;
  }
}

/***************************************************/
/* LEGENDRE transforms */
/* for Arpege and ECMWF global fields in spherical harmonics */
/***************************************************/
/* 1. initialise associated legendre functions for a single latitude */
/* 'legendre' is a vector of length (NSMAX+1)*(NSMAX+2)/2 
   What is the best way to treat a list of latitudes? As rows or columns?
   ATTENTION: it is impossible to hold really large sets in memory
   For HRES (IFS), e.g. NSMAX=o(1000) : ~500000 values/lat = 4MB/lat -> need >4GB!
   So then you have to calculate 1 lat at a time, on the fly.
*/
#define klf(n,m) ( ((n+1)*(n))/2 + m )
void legendre_init1(double *x, int *NSMAX,double* legendre){
  int n,m;
  double bn1,bn,s,mu;
  double SQ3=sqrt(3.00),SQ15=sqrt(1.5);

  mu = *x;

  s=sqrt(1- mu*mu);
  legendre[klf(0,0)] = 1.;
  legendre[klf(1,0)] = SQ3* mu;
  legendre[klf(1,1)] = SQ15* s;

  bn=SQ3;
  for(n=2;n<*NSMAX;n++){
    bn1=bn;
    bn=sqrt( (double) (4*n*n-1))/n;

    legendre[klf(n,0)] = (legendre[klf(n-1,0)]*mu-bn1*legendre[klf(n-2,0)])/bn;
    legendre[klf(n,1)] = sqrt(( (double) (2*n+1))/(2*n-1))*legendre[klf(n-1,0)] 
                         + sqrt(n/(n+1.))*legendre[klf(n,0)]*mu/s;
    legendre[klf(n,n)] = sqrt(((double) (2*n+1))/(2*n))*legendre[klf(n-1,n-1)]*s;

    for(m=2;m<=n;m++){
      legendre[klf(n,m)] = sqrt(((double) ((2*m-1)*(n+m-1)*(n+m-3)))/ ((2*n-3)*(n+m)*(n+m-2)) )*legendre[klf(n-2,m-2)]
                         - sqrt(((double) ((2*n+1)*(n+m-1)*(n-m+1))) /((2*n-1)*(n+m)*(n+m-2)) )*legendre[klf(n-1,m-2)]*mu
                         + sqrt(((double) ((2*n+1)*(n-m)))/((2*n-1)*(n+m)) )*legendre[klf(n-1,m)]*mu;
    }
  }
}

/* faster, different ordering ...*/
/* directly in the ordering of ECMWF grib files */
void legendre_init2(double *lat, int *NSMAX,double* legendre){
  int n,m,pos;
  double b1,b2,b3,b4,zn,zm,z1,z2,zsin,zcos;

  zsin=sin(*lat);
  zcos=sqrt(1.-zsin*zsin);

  legendre[0] = 1.;

  z1=sqrt(3.0);
  legendre[1] = z1 * zsin;

  pos=1;
  for(m=0;m< *NSMAX;m++){
    zm = (double)m;
    b1 = sqrt(2.*zm+3.);
    b2 = 1./sqrt(2.*zm);
        
    if(m){ 
      z2 = z1*zcos*b2;
      z1 = z2*b1;

      pos++;
      legendre[pos] = z2;
      pos++;
      legendre[pos] = z1 * zsin;
    }

    b4=1./b1; /* b4 for n=m+2: sqrt(1./(2m+3)) */
/* check: is n up to NSMAX+1???*/
    for(n=m+2;n<= *NSMAX;n++){
      zn=n;
      b3=sqrt((4.*zn*zn-1.)/(zn*zn-zm*zm));
      pos++;
      legendre[pos] = b3*(zsin*legendre[pos-1] - b4*legendre[pos-2]);
      b4=1./b3;
    }
  }
}
/* Arpege and ecmwf use different ordering. So we re-order Arpege data,. */
/* this code may be better paced in Rfa, but maybe there are other uses... */
void legendre_reorder(double *indata,double *outdata, int *nsmax){
  int n,m,pos,MAX;
  MAX=*nsmax+1;
  pos=0;
  for(m=0;m < MAX;m++){
    for(n=m;n<MAX;n++){
      outdata[pos++]=indata[(n*(n+1))/2+m];
    }
  }
}

/* converting back to Arpege ordering */
void legendre_reorder_inv(double *indata,double *outdata, int *nsmax){
  int n,m,pos,MAX;
  MAX=*nsmax+1;
  pos=0;
  for(m=0;m < MAX;m++){
    for(n=m;n<MAX;n++){
      outdata[(n*(n+1))/2+m] = indata[pos++];
    }
  }
}




