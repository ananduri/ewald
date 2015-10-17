#include "ewald.h"

double selfint(double alpha, double cellsize, int bsize){
	
	return -2*D*rnn*rnn*rnn*bsize*cellsize*cellsize*cellsize*alpha*alpha*alpha/sqrt(M_PI)/3;
}
