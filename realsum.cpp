#include "ewald.h"

double realsum(double x, double y, double z, int m, double alpha, int real_cut, double cellsize, int bsize, double *totNNenergy, double* intmat){

	extern inline double B(double, double);
	extern inline double C(double, double);
	
	double real=0;
	
	double off[12] = {0, 0, 0,
			0, 0.25, 0.25,
			0.25, 0, 0.25,
			0.25, 0.25, 0};

	double dipoff[12] = {is3, is3, is3,
				-is3, is3, is3,
				is3, -is3, is3,
				is3, is3, -is3};

	double* mu1;
	double* mu2;

	double* boffset1;

	double indenergy;
	double NNenergy;

	int N = bsize*cellsize*cellsize*cellsize;

	boffset1 = &off[3*m];

	for(int u=0; u<cellsize; ++u){
		for(int v=0; v<cellsize; ++v){
			for(int w=0; w<cellsize; ++w){
				for(int s=0; s<bsize; s++){

	double* boffset2;
	boffset2 = &off[3*s];

	double r, X, Y, Z;
	double first, second, dot;

	indenergy = 0;
	NNenergy = 0;

	for (int i=(-1)*real_cut; i<=real_cut; ++i){
		for (int j=(-1)*real_cut; j<=real_cut; ++j){
			for (int k=(-1)*real_cut; k<=real_cut; ++k){

				X = 0.5*((i+k)*cellsize + x-u + z-w) + boffset1[0] - boffset2[0]; 
				Y = 0.5*((i+j)*cellsize + x-u + y-v) + boffset1[1] - boffset2[1];
				Z = 0.5*((j+k)*cellsize + y-v + z-w) + boffset1[2] - boffset2[2];

				r = sqrt(X*X + Y*Y + Z*Z);

				if(r>0.001){

					mu1 = &dipoff[3*m];
	
					mu2 = &dipoff[3*s];
					
					dot = (mu1[0]*mu2[0] + mu1[1]*mu2[1] + mu1[2]*mu2[2]);

					if(r < 0.4) { 
					//if((r < 0.4) && (i==0) && (j==0) && (k==0)) { //which is right here? actually this is physically wrong?
						NNenergy += J*dot/2;
					}

					first = dot * B(r,alpha);

					second = -(mu1[0]*(X) + mu1[1]*(Y) + mu1[2]*(Z)) * (mu2[0]*(X) + mu2[1]*(Y) + mu2[2]*(Z)) * C(r,alpha);

					real += rnn*rnn*rnn*D*0.5*(first + second);
					
					indenergy += rnn*rnn*rnn*D*0.5*(first+second);
				}
			}
		}
	}

		intmat[(int)((bsize*cellsize*cellsize*x + bsize*cellsize*y + bsize*z + m)*N + bsize*cellsize*cellsize*u + bsize*cellsize*v + bsize*w + s)] += indenergy + NNenergy;
		*totNNenergy += NNenergy;
			}
			}
		}
	}
	
	return real;
}
