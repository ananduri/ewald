#include "ewald.h"

void dipassign(double* dip, int m, bool ori) {

	dip[0] = 1.0;
	dip[1] = 1.0;
	dip[2] = 1.0;

	if(m==1) {dip[0] *= -1;}
	else if(m==2) {dip[1] *= -1;}
	else if(m==3) {dip[2] *= -1;}
	
	//normalization and orientation
	for(int x=0; x<3; x++) {
		dip[x] = (ori) ? dip[x] : -dip[x]; //better would be dip[x] -= (!ori)*2*dip[x];
		dip[x] /= sqrt(3);
	}
	
	return;
}	
