#include "ewald.h"

double getenergy(double* intmat, int* state, int N) {
	
	double sum1;
	double sum2=0;

	for(int i =0; i<N; i++) {
		sum1 = 0;
		for(int j = 0; j<N; j++) {
			sum1 += intmat[i*N + j] * state[j];
		}

		sum2 += state[i] * sum1;
	}
	return sum2;
}
