#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;
#ifndef ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include "proj.h"
#else
#include "proj_api.h"
#endif
// [[Rcpp::export]]
NumericVector proj_version() {
  
  NumericVector result(3) ;
#ifndef ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
  result(0) = PROJ_VERSION_MAJOR ;
  result(1) = PROJ_VERSION_MINOR ;
  result(2) = PROJ_VERSION_PATCH ;
#else
  result(0) = floor(PJ_VERSION/100.) ;
  result(1) = floor( (PJ_VERSION % 100) / 10.) ;
  result(2) = PJ_VERSION % 10 ;
#endif
  return result;
}
