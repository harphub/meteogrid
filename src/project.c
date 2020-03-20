/* Simple wrapper script for proj[4] routines. */

#include <stdio.h>
#include "R.h"

/* new interface: TO DO */
#ifdef PROJ5
#include "proj.h"
void Rproj(double *u, double *v, int *npoints,
           char *parms,
           int *inverse){
  PJ *P;
  PJ_UV data;
  int i;

  if ( !(P=proj_create(XXX, parms)) ){
    Rprintf("ERROR: Projection Initialisation Error in\n");
    Rprintf("%s ", *parms);
    Rprintf("\n");
    Rf_error("Can not initialize projection.\n");
  }
  for(i=0; i < *npoints ; i++){
    if(R_FINITE(data.u = u[i]) && R_FINITE(data.v = v[i])) {
      data = (*inverse) ? pj_inv(data,ref) : pj_fwd(data,ref) ;
      if(data.u == HUGE_VAL) u[i] = v[i] = NA_REAL;
      else {
        u[i] = data.u;
        v[i] = data.v;
      }
    }
  }
  proj_destroy(P);
}

#endif

/* old interface (OK until 2020) */
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include "proj_api.h"

void Rproj4(double *u,double *v,int *npoints,char **parms,
            int *nparms,int *inverse){
  projPJ *ref;
  projUV data;
  int i;

  if ( !(ref=pj_init(*nparms,parms)) ){
    Rprintf("ERROR: Projection Initialisation Error in\n");
    for(i=0;i<*nparms;i++) Rprintf("%s ",parms[i]);
    Rprintf("\n");
    Rf_error("Can not initialize projection.\n");
  }
  for(i=0;i< *npoints;i++){
    if(R_FINITE(data.u = u[i]) && R_FINITE(data.v = v[i])) {
      data = (*inverse) ? pj_inv(data,ref) : pj_fwd(data,ref) ;
      if(data.u == HUGE_VAL) u[i] = v[i] = NA_REAL;
      else {
        u[i] = data.u;
        v[i] = data.v;
      }
    }
  }
  pj_free(ref);
}


