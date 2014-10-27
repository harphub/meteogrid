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

