#include "ewald.h"

double realsum(double x, double y, double z, double alpha, int recip_cut, int real_cut, double cellsize, bool charray[]){
	double real=0;

	double mu1[3];
	double mu2[3];

	for(int u=0; u<cellsize; ++u){
		for(int v=0; v<cellsize; ++v){
			for(int w=0; w<cellsize; ++w){

	//do sum in real-space
	double r, X, Y, Z;
	double first, second;
	for (int i=(-1)*real_cut; i<=real_cut; ++i){
		for (int j=(-1)*real_cut; j<=real_cut; ++j){
			for (int k=(-1)*real_cut; k<=real_cut; ++k){
				X = ((double) i)*cellsize + x-u;
				Y = ((double) j)*cellsize + y-v;
				Z = ((double) k)*cellsize + z-w;
				r = sqrt(X*X + Y*Y + Z*Z);

				if(r>0.001){

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


					first = (mu1[0]*mu2[0] + mu1[1]*mu2[1] + mu1[2]*mu2[2]) * B(r,alpha);

					second = -(mu1[0]*(X) + mu1[1]*(Y) + mu1[2]*(Z)) * (mu2[0]*(X) + mu2[1]*(Y) + mu2[2]*(Z)) * C(r,alpha);

					//std::cout << "C = " << C(r,alpha) << "\n";

					real += 0.5*(first + second);
					}
				}
			}
		}

			}
		}
	}
	return real;
}
