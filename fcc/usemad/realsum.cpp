#include "ewald.h"

double realsum(double x, double y, double z, int m, int p, double alpha, int real_cut, double cellsize, int bsize, double *totNNenergy, double* intmat, double* Jmat, double* Dmat){

	int fsize = 4;

	extern inline double B(double, double);
	extern inline double C(double, double);
	extern inline double A(double, double);
	
	double real=0;
	
	double off[12] = {0, 0, 0,
			0, 0.25, 0.25,
			0.25, 0, 0.25,
			0.25, 0.25, 0};


	//double off[12] = {0, 0, 0,
	//		0.5, 0, 0};

	double dipoff[12] = {is3, is3, is3,
				-is3, is3, is3,
				is3, -is3, is3,
				is3, is3, -is3};

	
	double foff[12] = {0, 0, 0,
			0, 0.5, 0.5,
			0.5, 0, 0.5,
			0.5, 0.5, 0};

	double* mu1;
	double* mu2;

	double* boffset1;
	double* foffset1;

	double indenergy;
	double NNenergy;

	int N = fsize*bsize*cellsize*cellsize*cellsize;

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

	double r, X, Y, Z;
	double first, second, dot;
	double ch1, ch2;

	indenergy = 0;
	NNenergy = 0;
	
	double Jdis,Ddis;

	for (int i=(-1)*real_cut; i<=real_cut; ++i){
		for (int j=(-1)*real_cut; j<=real_cut; ++j){
			for (int k=(-1)*real_cut; k<=real_cut; ++k){

				X = ((double) i)*cellsize + x-u + boffset1[0]-boffset2[0] + foffset1[0]-foffset2[0];	
				Y = ((double) j)*cellsize + y-v + boffset1[1]-boffset2[1] + foffset1[1]-foffset2[1];	
				Z = ((double) k)*cellsize + z-w + boffset1[2]-boffset2[2] + foffset1[2]-foffset2[2];	
		
				r = sqrt(X*X + Y*Y + Z*Z);

				if(r>0.001){

					mu1 = &dipoff[3*m];
	
					mu2 = &dipoff[3*s];
					
					dot = (mu1[0]*mu2[0] + mu1[1]*mu2[1] + mu1[2]*mu2[2]);

					if(r < 0.4) { 
						Jdis = Jmat[(int)((fsize*bsize*cellsize*cellsize*x + fsize*bsize*cellsize*y + fsize*bsize*z + fsize*m + p)*N + fsize*bsize*cellsize*cellsize*u + fsize*bsize*cellsize*v + fsize*bsize*w + fsize*s + q)];
						NNenergy += Jdis*dot/2;
					}

					first = dot * B(r,alpha);

					second = -(mu1[0]*(X) + mu1[1]*(Y) + mu1[2]*(Z)) * (mu2[0]*(X) + mu2[1]*(Y) + mu2[2]*(Z)) * C(r,alpha);

					Ddis = Dmat[(int)((fsize*bsize*cellsize*cellsize*x + fsize*bsize*cellsize*y + fsize*bsize*z + fsize*m + p)*N + fsize*bsize*cellsize*cellsize*u + fsize*bsize*cellsize*v + fsize*bsize*w + fsize*s + q)];

					real += rnn*rnn*rnn*Ddis*0.5*(first + second); //factor of half is because summing over whole matrix, instead of just upper triangle
					
					indenergy += rnn*rnn*rnn*Ddis*0.5*(first+second);

					//ch1 = 2*(m==1) - 1;
					//ch2 = 2*(s==1) - 1;
					
					//real += ch1*ch2*0.5*A(r,alpha); //including factor of half from whole matrix instead of upper half
					//indenergy += ch1*ch2*0.5*A(r,alpha);
				}
			}
		}
	}

		intmat[(int)((fsize*bsize*cellsize*cellsize*x + fsize*bsize*cellsize*y + fsize*bsize*z + fsize*m + p)*N + fsize*bsize*cellsize*cellsize*u + fsize*bsize*cellsize*v + fsize*bsize*w + fsize*s + q)] += indenergy + NNenergy; 
		//intmat[(int)((fsize*bsize*cellsize*cellsize*x + fsize*bsize*cellsize*y + fsize*bsize*z + fsize*m + p)*1 + N*(fsize*bsize*cellsize*cellsize*u + fsize*bsize*cellsize*v + fsize*bsize*w + fsize*s + q))] += indenergy + NNenergy; //had to reverse indexing

		*totNNenergy += NNenergy;
			}
			}
			}
		}
	}
	
	return real;
}
