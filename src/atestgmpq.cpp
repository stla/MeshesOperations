#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

#include "gmp.h"

std::string q2str(CGAL::Gmpq r){
	CGAL::Gmpz numer = r.numerator();
	CGAL::Gmpz denom = r.denominator();
	size_t n = numer.size();
	size_t d = denom.size();
	char *snumer = new char [n];
	char *sdenom = new char [d];
	mpz_get_str(snumer, 10, numer.mpz());
	mpz_get_str(sdenom, 10, denom.mpz());
	return snumer + "/" + sdenom;
}

// CGAL::Gmpq str2q(Rcpp::CharacterVector x){
// 	CGAL::Gmpq r(CGAL::Gmpq(Rcpp::as<std::string>(x)));
// 	return r;
// }

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
	std::string y = q2str(ggg):
	return y;
}

// [[Rcpp::export]]
double gtest3(Rcpp::CharacterVector x){
	CGAL::Gmpq gg(CGAL::Gmpq(Rcpp::as<std::string>(x(0))));
	double y = gg.to_double();
	return y;
}

