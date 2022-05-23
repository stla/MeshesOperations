#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

#include <gmp.h>

// [[Rcpp::export]]
double gtest(Rcpp::CharacterVector x){
	mpq_t gg;
	mpq_init_set_str(gg, "5/2", 10);
	double y = mpq_get_d(gg);
	return y;
}
