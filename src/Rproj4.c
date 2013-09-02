/* Simple wrapper script for proj4 routines. */
/* Can be optimised...                       */

#include <stdio.h>
#include "R.h"

#include "proj_api.h"


void Rproj4(double*,double*,int*,char**,int*,int*);


void Rproj4(double* u,double* v,int*npoints,char** parms,int* nparms,int* inverse){
  projPJ *ref;
  projUV data;
  int i;

  if ( !(ref=pj_init(*nparms,parms)) ){
    printf("Error: Projection Initialisation Error in\n");
    for(i=0;i<*nparms;i++) printf("%s ",parms[i]);
    exit(1);
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
}
