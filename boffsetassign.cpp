#include "ewald.h"
//use a switch
void boffsetassign(double* boff, int m) {
	
	if(m==0) {
		boff[0] = 0;
		boff[1] = 0;
		boff[2] = 0;
	}
	else if(m==1) {
		boff[0] = 0;
		boff[1] = 0.25;
		boff[2] = 0.25;
	}
	else if(m==2) { 
		boff[0] = 0.25;
		boff[1] = 0;
		boff[2] = 0.25;
	}
	else if(m==3) {
		boff[0] = 0.25;
		boff[1] = 0.25;
		boff[2] = 0;
	}

	return;
}
