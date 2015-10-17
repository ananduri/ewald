#include "ewald.h"

double recsum(double x, double y, double z, double m, double alpha, int recip_cut, double cellsize, int bsize, bool charray[], double* intmat){
	using namespace std;	

	double recip=0, k2;
	
	double indenergy;

	int N = bsize*cellsize*cellsize*cellsize;

	double mu1[3];
	double mu2[3];

	double dipk1;
	double dipk2;

	double boffset1[3] = {0.0}; //try this initialization with boffset2 as well
	boffsetassign(boffset1,m);

	/*if(m==0) { 
		boffset1[0] = 0;
		boffset1[1] = 0;
		boffset1[2] = 0;
	}
	else if(m==1) {
		boffset1[0] = 0.25;
		boffset1[1] = 0.25;
		boffset1[2] = 0;
	}
	else if(m==2) { 
		boffset1[0] = 0;
		boffset1[1] = 0.25;
		boffset1[2] = 0.25;
	}
	else if(m==3) { 
		boffset1[0] = 0.25;
		boffset1[1] = 0;
		boffset1[2] = 0.25;
	}*/


	for(int u=0; u<cellsize; ++u){
		for(int v=0; v<cellsize; ++v){
			for(int w=0; w<cellsize; ++w){
				for(int s=0; s<bsize; s++){

	double boffset2[3];
	boffsetassign(boffset2,s);

	/*if(s==0) { 
		boffset2[0] = 0;
		boffset2[1] = 0;
		boffset2[2] = 0;
	}
	else if(s==1) {
		boffset2[0] = 0.25;
		boffset2[1] = 0.25;
		boffset2[2] = 0;
	}
	else if(s==2) { 
		boffset2[0] = 0;
		boffset2[1] = 0.25;
		boffset2[2] = 0.25;
	}
	else if(s==3) { 
		boffset2[0] = 0.25;
		boffset2[1] = 0;
		boffset2[2] = 0.25;
	}*/


	double kterm;

	double I,J,K;

	double sfx,sfy,sfz;

	indenergy=0;

	//do sum in k-space (leaving out k=0 term)
	for (int i=(-1)*recip_cut; i<=recip_cut; i++){
		for (int j=(-1)*recip_cut; j<=recip_cut; j++){
			for (int k=(-1)*recip_cut; k<=recip_cut; k++){

				I = (double) (i - j + k); //FCC
				J = (double) (i + j - k);
				K = (double) (-i + j + k);

				sfx = 0.5*(x-u + z-w) + boffset1[0] - boffset2[0];
				sfy = 0.5*(x-u + y-v) + boffset1[1] - boffset2[1];
				sfz = 0.5*(y-v + z-w) + boffset1[2] - boffset2[2];


				k2 =  (I*I + J*J + K*K)/(cellsize*cellsize);
					
				if (k2 > 0.01/cellsize/cellsize){


					dipassign(mu1, m, charray[(int)(bsize*cellsize*cellsize*x + bsize*cellsize*y + bsize*z + m)]);

					dipassign(mu2, s, charray[(int)(bsize*cellsize*cellsize*u + bsize*cellsize*v + bsize*w + s)]);


					/*if(charray[(int)(bsize*cellsize*cellsize*x + bsize*cellsize*y + bsize*z + m)]){
						mu1[0] = 0;
						mu1[1] = 0;
						mu1[2] = 1;
					}
					else {
						mu1[0] = 0;
						mu1[1] = 0;
						mu1[2] = -1;
					}
					if(charray[(int)(bsize*cellsize*cellsize*u + bsize*cellsize*v + bsize*w + s)]){
						mu2[0] = 0;
						mu2[1] = 0;
						mu2[2] = 1;
					}
					else {
						mu2[0] = 0;
						mu2[1] = 0;
						mu2[2] = -1;
					}*/

					dipk1 = mu1[0]*I/cellsize + mu1[1]*J/cellsize + mu1[2]*K/cellsize;
					dipk2 = mu2[0]*I/cellsize + mu2[1]*J/cellsize + mu2[2]*K/cellsize;

					kterm = 1;
					kterm *= dipk1*dipk2*4*M_PI*M_PI;
					kterm *= 1.0/(M_PI*k2);
					kterm *= structure2(sfx, sfy, sfz, I/cellsize,J/cellsize,K/cellsize);
					kterm *= exp(-M_PI*M_PI*k2/(alpha*alpha));
					kterm /= cellsize*cellsize*cellsize;

					recip += rnn*rnn*rnn*D*4.0*0.5*kterm;

					indenergy += rnn*rnn*rnn*D*4.0*0.5*kterm;
					}
				}
			}
		}

		intmat[(int)((bsize*cellsize*cellsize*x + bsize*cellsize*y + bsize*z + m)*N + bsize*cellsize*cellsize*u + bsize*cellsize*v + bsize*w + s)] += indenergy;
		}
		}
		}
		}

	//cout << "recenergy: " << recip << '\n';

	
	return recip;
}
