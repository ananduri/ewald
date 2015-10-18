#include "ewald.h"

void foffsetassign(double* boff, int m) {
	
	double off[12] = {0, 0, 0,
			0, 0.5, 0.5,
			0.5, 0, 0.5,
			0.5, 0.5, 0};

	//should probably put some checks here
	boff = &off[3*m];

	/*if(m==0) {
		boff[0] = 0;
		boff[1] = 0;
		boff[2] = 0;
	}
	else if(m==1) {
		boff[0] = 0;
		boff[1] = 0.5;
		boff[2] = 0.5;
	}
	else if(m==2) { 
		boff[0] = 0.5;
		boff[1] = 0;
		boff[2] = 0.5;
	}
	else if(m==3) {
		boff[0] = 0.5;
		boff[1] = 0.5;
		boff[2] = 0;
	}*/

	return;
}
