#include "ewald.h"

double recsum(double x, double y, double z, int m, int p, double alpha, int recip_cut, double cellsize, int bsize, double* intmat, double* Dmat){

	int fsize = 4;

	extern inline double structure2(double, double, double, double, double, double);

	double recip=0, k2;
	
	double indenergy;

	int N = fsize*bsize*cellsize*cellsize*cellsize;

	double* mu1;
	double* mu2;

	double dipk1;
	double dipk2;

	//double off[12] = {0, 0, 0,
	//		0.5, 0, 0};

	double off[12] = {0, 0, 0,
			0, 0.25, 0.25,
			0.25, 0, 0.25,
			0.25, 0.25, 0};

	double foff[12] = {0, 0, 0,
			0, 0.5, 0.5,
			0.5, 0, 0.5,
			0.5, 0.5, 0};

	double dipoff[12] = {is3, is3, is3,
				-is3, is3, is3,
				is3, -is3, is3,
				is3, is3, -is3};

	/*double fkdict[12] = {0, 0, 0,
				1, 0, 0,
				0, 1, 0,
				0, 0, 1}; */

	double* boffset1;
	double* foffset1;
	boffset1 = &off[3*m];
	foffset1 = &foff[3*p];

	for(int u=0; u<cellsize; ++u){ 
		for(int v=0; v<cellsize; ++v){
			for(int w=0; w<cellsize; ++w){
				for(int s=0; s<bsize; s++){
				for(int q=0; q<fsize; q++){

	double* boffset2;
	double* foffset2;
	boffset2 = &off[3*s];
	foffset2 = &foff[3*q];

	double kterm;
	double mkterm;

	double I,J,K;

	double sfx,sfy,sfz;
	double ch1,ch2;
	
	double* fkvec;

	indenergy=0;

	double Jdis,Ddis;

	//do sum in k-space (leaving out k=0 term) //make sure cellsize is even
	for (int i=(-1)*recip_cut; i<=recip_cut; i++){
		for (int j=(-1)*recip_cut; j<=recip_cut; j++){
			for (int k=(-1)*recip_cut; k<=recip_cut; k++){

				//factors of 1/2, 2pis, and cellsize were taken from here and put in below, so be careful
				/*I = (double) (i - j + k); //FCC
				J = (double) (i + j - k);
				K = (double) (-i + j + k);*/

				//here define the allowed kvecs
				//currently: factor of 2pi/cellsize included LATER
				I = (double)i;// + cellsize*fkvec[0];
				J = (double)j;// + cellsize*fkvec[1];
				K = (double)k;// + cellsize*fkvec[2]; 

				//sfx = 0.5*(x-u + z-w) + boffset1[0] - boffset2[0];
				//sfy = 0.5*(x-u + y-v) + boffset1[1] - boffset2[1];
				//sfz = 0.5*(y-v + z-w) + boffset1[2] - boffset2[2];

				sfx = x-u + boffset1[0]-boffset2[0] + foffset1[0]-foffset2[0];
				sfy = y-v + boffset1[1]-boffset2[1] + foffset1[1]-foffset2[1];
				sfz = z-w + boffset1[2]-boffset2[2] + foffset1[2]-foffset2[2];

				k2 =  (I*I + J*J + K*K)/(cellsize*cellsize); //still not including the (2pi)^2
					
				if (k2 > 0.01/cellsize/cellsize){

					mu1 = &dipoff[3*m];

					mu2 = &dipoff[3*s];

					//here doing a dot product, of mu with k, including cellsize but not 2pi in k
					dipk1 = mu1[0]*I/cellsize + mu1[1]*J/cellsize + mu1[2]*K/cellsize;
					dipk2 = mu2[0]*I/cellsize + mu2[1]*J/cellsize + mu2[2]*K/cellsize;

					/*kterm = 1;
					kterm *= dipk1*dipk2*4*M_PI*M_PI;
					kterm *= 1.0/(M_PI*k2);
					kterm *= structure2(sfx, sfy, sfz, I/cellsize,J/cellsize,K/cellsize);
					kterm *= exp(-M_PI*M_PI*k2/(alpha*alpha));
					kterm /= cellsize*cellsize*cellsize;

					kterm /= 16.0; //where the fuck does this come from, possibly a factor of 4 from fkvec duplicating?
					//kterm /= 4; */

					kterm = 1;
					kterm *= dipk1*dipk2;
					kterm /= 2*M_PI*k2*cellsize*cellsize*cellsize;
					kterm *= structure2(sfx,sfy,sfz, I/cellsize, J/cellsize, K/cellsize);
					kterm *= exp(-M_PI*M_PI*k2/(alpha*alpha));

					kterm *= 4*M_PI*M_PI; //also in dipole dot prods

					Ddis = Dmat[(int)((fsize*bsize*cellsize*cellsize*x + fsize*bsize*cellsize*y + fsize*bsize*z + fsize*m + p)*N + fsize*bsize*cellsize*cellsize*u + fsize*bsize*cellsize*v + fsize*bsize*w + fsize*s + q)];

					recip += rnn*rnn*rnn*Ddis*kterm; 

					indenergy += rnn*rnn*rnn*Ddis*kterm;
					
					/*ch1 = 2*(m==1) - 1;
					ch2 = 2*(s==1) - 1;
	
					mkterm = 1;
					mkterm *= ch1*ch2;
					mkterm /= 2*M_PI*k2*cellsize*cellsize*cellsize; 
					mkterm *= structure2(sfx,sfy,sfz, I/cellsize, J/cellsize, K/cellsize);
					mkterm *= exp(-M_PI*M_PI*k2/(alpha*alpha));
	
					recip += mkterm;
					indenergy += mkterm; //including 0.5 because doing whole matrix instead of half (even with no 0.5 written)*/
				}
				}
			}
		}

		intmat[(int)((fsize*bsize*cellsize*cellsize*x + fsize*bsize*cellsize*y + fsize*bsize*z + fsize*m + p)*N + fsize*bsize*cellsize*cellsize*u + fsize*bsize*cellsize*v + fsize*bsize*w + fsize*s + q)] += indenergy;
		//intmat[(int)((fsize*bsize*cellsize*cellsize*x + fsize*bsize*cellsize*y + fsize*bsize*z + fsize*m + p)*1 + N*(fsize*bsize*cellsize*cellsize*u + fsize*bsize*cellsize*v + fsize*bsize*w + fsize*s + q))] += indenergy; //had to reverse indices
		}
		}
		}
		}
		}
	return recip;
}
