#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

#include <gmpxx.h>

// [[Rcpp::export]]
double gtest(Rcpp::CharacterVector x){
	mpq_t gg;
	mpq_init_set_str(gg, x(0), 10);
	double y = mpq_get_d(gg);
	return y;
}
