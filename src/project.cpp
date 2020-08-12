#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

#ifdef PROJ5
#include "proj.h"
#else
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include "proj_api.h"
#endif

// [[Rcpp::export]]
NumericMatrix mg_project( NumericVector x, NumericVector y, std::string proj_string,
            bool inverse=false){

#ifdef PROJ5
  PJ *P;
  PJ_COORD data;
#else
  projPJ P;
  projUV data;
#endif

  int i, npoints = x.length();
  const char *parms = proj_string.c_str();
  NumericMatrix result(npoints, 2);

#ifdef PROJ5
  P = proj_create(PJ_DEFAULT_CTX, *parms);
#else
  P = pj_init_plus(parms);
#endif

  if ( !P ) {
    Rprintf("ERROR: Projection Initialisation Error in\n");
    Rprintf("%s\n", *parms);
    Rf_error("Can not initialize projection.\n");
  }

#ifdef PROJ5
  if (! inverse) {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(data.lp.lam = proj_torad(x[i])) &&
          R_FINITE(data.lp.phi = proj_torad(y[i]))) {
        data = proj_trans(P, PJ_FWD, data) ;
        if (data.lp.lam == HUGE_VAL) result[i, 0] = result[i, 1] = NA_REAL;
        else {
          result[i, 0] = data.xy.x;
          result[i, 1] = data.xy.y;
        }
      }
    }
  }
  else {
    for (i=0; i < *npoints ; i++) {
      if (R_FINITE(data.xy.x = x[i]) && R_FINITE(data.xy.y = y[i])) {
        data = proj_trans(P, PJ_INV, data) ;
        if (data.lp.lam == HUGE_VAL) result[i, 0] = result[i, 1] = NA_REAL;
        else {
          result[i, 0] = proj_todeg(data.lp.lam);
          result[i, 1] = proj_todeg(data.lp.phi);
        }
      }
    }
  }
  proj_destroy(P);
#else
  if (! inverse) {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(data.u = DEG_TO_RAD * x[i]) &&
          R_FINITE(data.v = DEG_TO_RAD * y[i])) {
        data = pj_fwd(data, P) ;
        if (data.u == HUGE_VAL) result[i, 0] = result[i, 1] = NA_REAL;
        else {
          result[i, 0] = data.u;
          result[i, 1] = data.v;
        }
      }
    }
  }
  else {
    for (i=0; i < npoints ; i++) {
      if (R_FINITE(data.u = x[i]) && R_FINITE(data.v = y[i])) {
        data = pj_inv(data, P) ;
        if (data.u == HUGE_VAL) result[i, 0] = result[i, 1] = NA_REAL;
        else {
          result[i, 0] = data.u * RAD_TO_DEG;
          result[i, 1] = data.v * RAD_TO_DEG;
        }
      }
    }
  }
  pj_free(P);
#endif
  return result;
}

