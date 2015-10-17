#include "ewald.h"

double realsum(double x, double y, double z, double m, double alpha, int real_cut, double cellsize, int bsize, bool charray[], double *totNNenergy, double* intmat){
	using namespace std;
	
	double real=0;
	

	//double checkintmat = 0;

	double mu1[3];
	double mu2[3];

	double boffset1[3] = {0.0};

	double indenergy;
	double NNenergy;

	int N = bsize*cellsize*cellsize*cellsize;

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

	boffsetassign(boffset1, m);


	for(int u=0; u<cellsize; ++u){
		for(int v=0; v<cellsize; ++v){
			for(int w=0; w<cellsize; ++w){
				for(int s=0; s<bsize; s++){


	double boffset2[3];
	boffsetassign(boffset2, s);

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


					
					dot = (mu1[0]*mu2[0] + mu1[1]*mu2[1] + mu1[2]*mu2[2]);

					//if(r<1.0) {printf("r: %.3f, i: %d, j: %d, k: %d\n",r,i,j,k);}
				
					if((r < 0.4) && (i==0) && (j==0) && (k==0)) { //this and the below line are both OK now
					//if(r<0.4) 
						NNenergy += J*dot/2;
						//std::cout<< "NNplus: " << J*dot/2 << "\n";
					}


					first = dot * B(r,alpha);

					second = -(mu1[0]*(X) + mu1[1]*(Y) + mu1[2]*(Z)) * (mu2[0]*(X) + mu2[1]*(Y) + mu2[2]*(Z)) * C(r,alpha);

					real += rnn*rnn*rnn*D*0.5*(first + second);
					
					indenergy += rnn*rnn*rnn*D*0.5*(first+second);
					}
			}
		}
	}

		//std::cout << "NNenergy: " << NNenergy << "\n";

		intmat[(int)((bsize*cellsize*cellsize*x + bsize*cellsize*y + bsize*z + m)*N + bsize*cellsize*cellsize*u + bsize*cellsize*v + bsize*w + s)] += indenergy + NNenergy;
		*totNNenergy += NNenergy;
		//checkintmat += indenergy;

			}
			}
		}
	}
	
	//cout << "realenergy: " << real << '\n';
	//cout << "indenergy: " << checkintmat << '\n';
	
	return real;
}
