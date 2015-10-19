#include "ewald.h"

double selfint(double alpha, double cellsize){
	
	return -2*cellsize*cellsize*cellsize*alpha*alpha*alpha/sqrt(M_PI)/3;
}
