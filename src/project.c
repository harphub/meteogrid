/* Simple wrapper script for proj[4] routines. */

#include <stdio.h>
#include "R.h"

/* new interface: needs more validation */
//#define PROJ5
#ifdef PROJ5
#include "proj.h"
void Rproj(double *u, double *v, int *npoints,
           char **parms,
           int *inverse){
  PJ *P;
  PJ_COORD data;
  int i;

  if ( !(P=proj_create(PJ_DEFAULT_CTX, *parms)) ){
    Rprintf("ERROR: Projection Initialisation Error in\n");
    Rprintf("%s\n", *parms);
    Rf_error("Can not initialize projection.\n");
  }
  if (! *inverse) {
    for (i=0; i < *npoints ; i++) {
      if (R_FINITE(data.lp.lam = proj_torad(u[i])) &&
          R_FINITE(data.lp.phi = proj_torad(v[i]))) {
        data = proj_trans(P, PJ_FWD, data) ;
        if (data.lp.lam == HUGE_VAL) u[i] = v[i] = NA_REAL;
        else {
          u[i] = data.xy.x;
          v[i] = data.xy.y;
        }
      }
    }
  }
  else {
    for (i=0; i < *npoints ; i++) {
      if (R_FINITE(data.xy.x = u[i]) && R_FINITE(data.xy.y = v[i])) {
        data = proj_trans(P, PJ_INV, data) ;
        if (data.lp.lam == HUGE_VAL) u[i] = v[i] = NA_REAL;
        else {
          u[i] = proj_todeg(data.lp.lam);
          v[i] = proj_todeg(data.lp.phi);
        }
      }
    }
  }

  proj_destroy(P);
}

#else
/* old interface (OK until 2020) */
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include "proj_api.h"

void Rproj(double *u,double *v,int *npoints,char **parms,
            int *inverse){
  projPJ *P;
  projUV data;
  int i;

//  if ( !(ref=pj_init(*nparms,parms)) ){
  if ( !(P=pj_init_plus(*parms)) ){
    Rprintf("ERROR: Projection Initialisation Error in\n");
    Rprintf("%s\n", *parms);
    Rf_error("Can not initialize projection.\n");
  }
  if (! *inverse) {
    for (i=0; i < *npoints ; i++) {
      if (R_FINITE(data.u = DEG_TO_RAD*u[i]) &&
          R_FINITE(data.v = DEG_TO_RAD*v[i])) {
        data = pj_fwd(data, P) ;
        if (data.u == HUGE_VAL) u[i] = v[i] = NA_REAL;
        else {
          u[i] = data.u;
          v[i] = data.v;
        }
      }
    }
  }
  else {
    for (i=0; i < *npoints ; i++) {
      if (R_FINITE(data.u = u[i]) && R_FINITE(data.v = v[i])) {
        data = pj_inv(data, P) ;
        if (data.u == HUGE_VAL) u[i] = v[i] = NA_REAL;
        else {
          u[i] = data.u * RAD_TO_DEG;
          v[i] = data.v * RAD_TO_DEG;
        }
      }
    }
  }

  pj_free(P);
}

#endif

