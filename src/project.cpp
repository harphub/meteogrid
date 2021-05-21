#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;
// The "proj_api.h" version is deprecated, but still available up to PROJ.7 (2020)
// So the new API version is here for testing.
// Only people using PROJ.8 or later will really need the new interface.
// but you might as well use proj.h if you have v5.0 or later
#ifndef METEOGRID_OLD_PROJ_API
#include "proj.h"
#else
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include "proj_api.h"
#endif

// [[Rcpp::export]]
Rcpp::DataFrame mg_project( NumericVector x, NumericVector y, std::string proj_string,
            bool inverse=false){

  int i, npoints = x.length();
  const char *parms = proj_string.c_str();
  // NOTE: you could do inline replacements (i.e. use x and y for output, too)
  NumericVector result_x(npoints), result_y(npoints);

#ifndef METEOGRID_OLD_PROJ_API
  // The new interface (PROJ >= 5.0)
  PJ *P;
  PJ_COORD data;
  P = proj_create(PJ_DEFAULT_CTX, parms);

  if ( !P ) {
    Rprintf("ERROR: Projection Initialisation Error in\n");
    Rprintf("%s\n", *parms);
    Rf_error("Can not initialize projection.\n");
  }

  // NOTE: proj_torad(NA) seems to give 0
  if (! inverse) {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(x[i]) && R_FINITE(y[i])) {
        data.lp.lam = proj_torad(x[i]) ;
        data.lp.phi = proj_torad(y[i]) ;
        data = proj_trans(P, PJ_FWD, data) ;
        if (data.xy.x == HUGE_VAL || data.xy.y == HUGE_VAL) {
          result_x[i] = result_y[i] = NA_REAL;
        }
        else {
          result_x[i] = data.xy.x;
          result_y[i] = data.xy.y;
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  else {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(x[i]) && R_FINITE(y[i])) {
        data.xy.x = x[i] ;
        data.xy.y = y[i] ;
        data = proj_trans(P, PJ_INV, data) ;
        if (data.lp.lam == HUGE_VAL || data.lp.phi == HUGE_VAL) {
          result_x[i] = result_y[i] = NA_REAL;
        }
        else {
          result_x[i] = proj_todeg(data.lp.lam);
          result_y[i] = proj_todeg(data.lp.phi);
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  proj_destroy(P);

#else
  // old PROJ4 interface (proj_api,h)
  projPJ P;
  projUV data;
  P = pj_init_plus(parms);

  if ( !P ) {
    Rprintf("ERROR: Projection Initialisation Error in\n");
    Rprintf("%s\n", *parms);
    Rf_error("Can not initialize projection.\n");
  }
  if (! inverse) {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(x[i]) && R_FINITE(y[i])) {
        data.u = DEG_TO_RAD * x[i] ;
        data.v = DEG_TO_RAD * y[i] ;
        data = pj_fwd(data, P) ;
        if (data.u == HUGE_VAL || data.v == HUGE_VAL) {
          result_x[i] = result_y[i] = NA_REAL;
        }
        else {
          result_x[i] = data.u;
          result_y[i] = data.v;
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  else {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(x[i]) && R_FINITE(y[i])) {
        data.u = x[i] ;
        data.v = y[i] ;
        data = pj_inv(data, P) ;
        if (data.u == HUGE_VAL || data.v == HUGE_VAL) {
          result_x[i] = result_y[i] = NA_REAL;
        }
        else {
          result_x[i] = data.u * RAD_TO_DEG;
          result_y[i] = data.v * RAD_TO_DEG;
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  pj_free(P);
#endif

  return Rcpp::DataFrame::create(Named("x") = result_x, Named("y") = result_y);
}

