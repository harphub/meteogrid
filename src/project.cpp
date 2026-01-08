#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;
// The "proj_api.h" version is deprecated, but still available up to PROJ.7 (2020)
// So the new API version is here for testing.
// Only people using PROJ.8 or later will really need the new interface.
// but you might as well use proj.h if you have v5.0 or later
#ifndef ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include "proj.h"
#else
#include "proj_api.h"
#endif

// [[Rcpp::export]]
Rcpp::DataFrame mg_project( NumericVector x, NumericVector y, std::string proj_string,
            bool inverse=false){

  int i, npoints = x.length();
  const char *parms = proj_string.c_str();
  // NOTE: you could do inline replacements (i.e. use x and y for output, too)
  NumericVector result_x(npoints), result_y(npoints);

#ifndef ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
  // The new interface (PROJ >= 5.0)
  PJ *P;
  PJ_COORD input, output;
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
        input.lp.lam = proj_torad(x[i]) ;
        input.lp.phi = proj_torad(y[i]) ;
        output = proj_trans(P, PJ_FWD, input) ;
        if (output.xy.x == HUGE_VAL || output.xy.y == HUGE_VAL || proj_errno(P) > 0) {
          //Rprintf("ERROR %i %lf %lf\n", i, x[i], y[i]);
          //Rprintf("Resetting P: %i\n", proj_errno(P));
          proj_errno_reset(P);
          result_x[i] = result_y[i] = NA_REAL;
        }
        else {
          result_x[i] = output.xy.x;
          result_y[i] = output.xy.y;
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  else {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(x[i]) && R_FINITE(y[i])) {
        input.xy.x = x[i] ;
        input.xy.y = y[i] ;
        output = proj_trans(P, PJ_INV, input) ;
        if (output.lp.lam == HUGE_VAL || output.lp.phi == HUGE_VAL || proj_errno(P) > 0) {
          result_x[i] = result_y[i] = NA_REAL;
          proj_errno_reset(P);
        }
        else {
          result_x[i] = proj_todeg(output.lp.lam);
          result_y[i] = proj_todeg(output.lp.phi);
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  proj_destroy(P);

#else
  // old PROJ4 interface (proj_api,h)
  projPJ P;
  projUV input, output;
  P = pj_init_plus(parms);

  if ( !P ) {
    Rprintf("ERROR: Projection Initialisation Error in\n");
    Rprintf("%s\n", *parms);
    Rf_error("Can not initialize projection.\n");
  }
  if (! inverse) {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(x[i]) && R_FINITE(y[i])) {
        input.u = DEG_TO_RAD * x[i] ;
        input.v = DEG_TO_RAD * y[i] ;
        output = pj_fwd(input, P) ;
        if (output.u == HUGE_VAL || output.v == HUGE_VAL) {
          result_x[i] = result_y[i] = NA_REAL;
        }
        else {
          result_x[i] = output.u;
          result_y[i] = output.v;
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  else {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(x[i]) && R_FINITE(y[i])) {
        input.u = x[i] ;
        input.v = y[i] ;
        output = pj_inv(input, P) ;
        if (output.u == HUGE_VAL || output.v == HUGE_VAL) {
          result_x[i] = result_y[i] = NA_REAL;
        }
        else {
          result_x[i] = output.u * RAD_TO_DEG;
          result_y[i] = output.v * RAD_TO_DEG;
        }
      }
      else result_x[i] = result_y[i] = NA_REAL;
    }
  }
  pj_free(P);
#endif

  return Rcpp::DataFrame::create(Named("x") = result_x, Named("y") = result_y);
}

