#include <Rcpp.h>
using namespace Rcpp;
// from high resolution (or random points) to lower resolution
// by simple mean over all values within a grid box

// INPUT:
// px, py: vectors of length npoints (in grid coordinates!). Integer values
// pval: the original values
// gnx, gny : the dimensions of the (upscaling) target grid
// OUTPUT:
// gcount, gval: the number of points falling in every cell, and the mean value


// [[Rcpp::export]]
NumericMatrix upscale_by_mean(NumericVector px, NumericVector py, NumericVector pval,
                              int gnx, int gny)
{
  int i, j, npoints=px.length();
  NumericMatrix gval(gnx, gny), gcount(gnx, gny);
  // TODO: check that this is initialised to 0!
  //       check that gnx,gny are correctly set to integer

  //  Rprintf("Calculating.\n");
  for (i=0; i < npoints; i++) {
    // TODO: check for NA_REAL or NA_INTEGER ??
    if (px[i] != NA_REAL && py[i]!= NA_REAL && R_FINITE(pval[i]) ) {
      if ( px[i] > 0 && px[i] <= gnx && 
           py[i] > 0 && py[i] <= gny  ) {
        gcount[px[i]-1, py[i]-1] += 1;
        gval[px[i]-1, py[i]-1] += pval[i];
      }
    }
  }

//  Rprintf("Normalising.\n");
  for (i=0 ; i < gnx ; i++) 
    for (j=0 ; j < gny; j++) {
    if (gcount[i,j] > 0 ) gval[i,j] /= gcount[i,j];
    else gval[i,j] = NA_REAL;
  }
  return gval;
}

