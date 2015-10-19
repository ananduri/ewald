#include "ewald.h"

double B(double r, double alpha){
	double term1, term2;	
	
	term1 = (1 - erf(alpha*r))/(r*r*r);
	term2 = (2*alpha/sqrt(M_PI))*exp(-alpha*alpha*r*r)/(r*r);
	return term1+term2;
}

double C(double r, double alpha){
	double term1, term2;

	term1 = 3*(1 - erf(alpha*r))/(r*r*r*r*r);
	term2 = (2*alpha/sqrt(M_PI))*(2*alpha*alpha + 3/(r*r))*exp(-alpha*alpha*r*r)/(r*r);
	return term1+term2;
}