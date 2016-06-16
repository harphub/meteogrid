/*
#--------------------------------------#
# Part of R-package geogrid            #
# © Alex Deckmyn                       #
# RMI Belgium, alex.deckmyn@meteo.be   #
# Released under GPL-3 license         #
#--------------------------------------#
*/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include "R.h"
// from high resolution (or random points) to lower resolution
// by simple mean over all values within a grid box

void upscale_by_mean(int *npoints, int *px, int *py, double *pval,
                     int *gnx, int *gpx, int *gcount, double *gval);
// INPUT:
// npoints: number of points to be aggregated to grid
// px, py: vectors of length npoints (in grid coordinates!)
// pval: the original values
// gnx, gny : the dimensions of the (upscaling) target grid
// OUTPUT:
// gcount, gval: the number of points falling in every cell, and the mean value

void upscale_by_mean(int *npoints, int *px, int *py, double *pval,
                     int *gnx, int *gny, int *gcount, double *gval)
{
  int i;

//  Rprintf("npoints=%d newgrid: %d x %d \n",*npoints, *gnx, *gny);
//  for(i=0;i<10;i++) Rprintf("p%d: %d %d %lf\n", i+1, px[i], py[i], pval[i]);

//  Rprintf("Initialising.\n");
  for (i=0 ; i< *gnx * *gny; i++) gcount[i] = gval[i] = 0;
  
//  Rprintf("Calculating.\n");
  for (i=0; i<*npoints; i++) {
    if (R_FINITE(px[i]) && R_FINITE(py[i]) && R_FINITE(pval[i]) ) {
      if ( px[i] > 0 && px[i] <= *gnx && 
           py[i] > 0 && py[i] <= *gny  ) {
        gcount[ (py[i]-1) * *gnx + px[i] - 1] += 1;
        gval[(py[i]-1) * *gnx + px[i] - 1] += pval[i];
      }
    }
  }

//  Rprintf("Normalising.\n");
  for (i=0 ; i< *gnx * *gny; i++) {
    if (gcount[i] > 0 ) gval[i] /= gcount[i];
    else gval[i] = NA_REAL;
  }
}

