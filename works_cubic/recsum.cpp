#include "ewald.h"

double recsum(double x, double y, double z, double alpha, int recip_cut, int real_cut, double cellsize, bool charray[]){
	double recip=0, k2;

	double mu1[3];
	double mu2[3];

	double dipk1;
	double dipk2;

	
	for(int u=0; u<cellsize; ++u){
		for(int v=0; v<cellsize; ++v){
			for(int w=0; w<cellsize; ++w){


	double kterm;

	double I,J,K;

	//do sum in k-space (leaving out k=0 term)
	for (int i=(-1)*recip_cut; i<=recip_cut; i++){
		for (int j=(-1)*recip_cut; j<=recip_cut; j++){
			for (int k=(-1)*recip_cut; k<=recip_cut; k++){

				I = (double) i;
				J = (double) j;
				K = (double) k;

				k2 =  (I*I + J*J + K*K)/(cellsize*cellsize);
					
				if (k2 > 0.5/cellsize/cellsize){

					if(charray[(int)(cellsize*cellsize*x + cellsize*y + z)]){
						mu1[0] = 0;
						mu1[1] = 0;
						mu1[2] = 1;
					}
					else {
						mu1[0] = 0;
						mu1[1] = 0;
						mu1[2] = -1;
					}
					if(charray[(int)(cellsize*cellsize*u + cellsize*v + w)]){
						mu2[0] = 0;
						mu2[1] = 0;
						mu2[2] = 1;
					}
					else {
						mu2[0] = 0;
						mu2[1] = 0;
						mu2[2] = -1;
					}

					dipk1 = mu1[0]*I/cellsize + mu1[1]*J/cellsize + mu1[2]*K/cellsize;
					dipk2 = mu2[0]*I/cellsize + mu2[1]*J/cellsize + mu2[2]*K/cellsize;

					kterm = 1;
					kterm *= dipk1*dipk2*4*M_PI*M_PI;
					kterm *= 0.5/(M_PI*k2);
					kterm *= structure2(x-u,y-v,z-w,I/cellsize,J/cellsize,K/cellsize);
					kterm *= exp(-M_PI*M_PI*k2/(alpha*alpha));
					kterm /= cellsize*cellsize*cellsize;

					recip += kterm;
					}
				}
			}
		}

		}
		}
		}
	return recip;
}
