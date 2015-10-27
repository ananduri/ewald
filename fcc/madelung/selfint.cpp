#include "ewald.h"

double selfint(double alpha, double cellsize, int bsize){
	//fsize included below	
	double fsize = 4;
	//return -2*D*rnn*rnn*rnn*bsize*fsize*cellsize*cellsize*cellsize*alpha*alpha*alpha/sqrt(M_PI)/3;
	return -D*rnn*rnn*rnn*bsize*fsize*cellsize*cellsize*cellsize*alpha/sqrt(M_PI);
}
