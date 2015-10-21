#include "ewald.h"

double selfint(double alpha, double cellsize, int bsize){
	//fsize included below	
	return -2*D*rnn*rnn*rnn*bsize*bsize*cellsize*cellsize*cellsize*alpha*alpha*alpha/sqrt(M_PI)/3;
}
