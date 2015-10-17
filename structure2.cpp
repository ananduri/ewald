#include "ewald.h"

double structure2(double x, double y, double z, double i, double j, double k){
	return cos(M_PI*2.0*(x*((double) i) + y*((double) j) + z*((double) k)));
	}