#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

#include <gmpxx.h>

// [[Rcpp::export]]
double gtest(Rcpp::CharacterVector x){
	mpq_class gg == x(0);
	double y = mpq_class::mpq_get_d(gg);
	return y;
}
