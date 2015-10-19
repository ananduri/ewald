#include "ewald.h"

double recsum(double x, double y, double z, int m, double alpha, int recip_cut, double cellsize, int bsize, double* intmat){

	extern inline double structure2(double, double, double, double, double, double);

	using namespace std;	

	double recip=0, k2;
	
	double indenergy;

	int N = bsize*cellsize*cellsize*cellsize;

	double* mu1;
	double* mu2;

	double dipk1;
	double dipk2;

	double off[12] = {0, 0, 0,
			0, 0.25, 0.25,
			0.25, 0, 0.25,
			0.25, 0.25, 0};

	double dipoff[12] = {is3, is3, is3,
				-is3, is3, is3,
				is3, -is3, is3,
				is3, is3, -is3};

	double* boffset1;
	boffset1 = &off[3*m];

	for(int u=0; u<cellsize; ++u){
		for(int v=0; v<cellsize; ++v){
			for(int w=0; w<cellsize; ++w){
				for(int s=0; s<bsize; s++){

	double* boffset2;
	boffset2 = &off[3*s];

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

					mu1 = &dipoff[3*m];

					mu2 = &dipoff[3*s];

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
	
	return recip;
}
