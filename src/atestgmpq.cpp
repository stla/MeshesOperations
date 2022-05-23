#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif



// [[Rcpp::export]]
double gtest(Rcpp::CharacterVector x){
	CGAL::Gmpq gg(CGAL::Gmpq("5/2"));
	double y = gg.to_double();
	return y;
}

// [[Rcpp::export]]
std::string gtest2(Rcpp::CharacterVector x){
	CGAL::Gmpq gg(CGAL::Gmpq("5/2"));
	CGAL::Gmpq ggg = gg + 2;
	std::string y = ggg.numerator() + "/" + ggg.denominator();
	return y;
}

// [[Rcpp::export]]
double gtest3(Rcpp::CharacterVector x){
	CGAL::Gmpq gg(CGAL::Gmpq(Rcpp::as<std::string>(x(0))));
	double y = gg.to_double();
	return y;
}

