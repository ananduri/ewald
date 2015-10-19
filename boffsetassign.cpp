#include "ewald.h"

void boffsetassign(double* boff, int m) {


	double off[12] = {0, 0, 0,
			0, 0.25, 0.25,
			0.25, 0, 0.25,
			0.25, 0.25, 0}; //its because this guy is a local variable isnt it? when this function exits off is cleared from the stack and boff2 doesnt point to anything anymore. can test this

	//boff = &off[3*m];

	
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
