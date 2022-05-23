#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

#include <gmpxx.h>

// [[Rcpp::export]]
double gtest(Rcpp::CharacterVector x){
	mpq_class gg;
	gg = "5/2";
	double y = gg.get_d();
	return y;
}

// [[Rcpp::export]]
std::string gtest2(Rcpp::CharacterVector x){
	mpq_class gg;
	gg = "5/2";
	mpq_class ggg = gg + 2;
	std::string y = ggg.get_str();
	return y;
}

// [[Rcpp::export]]
double gtest3(Rcpp::CharacterVector x){
	mpq_class gg;
	gg = x(0);
	double y = gg.get_d();
	return y;
}

// [[Rcpp::export]]
double gtest4(Rcpp::CharacterVector x){
	mpq_class gg(x(0));
	double y = gg.get_d();
	return y;
}
