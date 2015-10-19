#include "ewald.h"

using namespace std;

int main(int argc, char *argv[]){

	double alpha;
	int real_cut, recip_cut, cellsize;
	
	cellsize = atoi(argv[1]);
	alpha = strtod(argv[2],NULL);
	real_cut = atoi(argv[3]);
	recip_cut = atoi(argv[4]);
	
	//for for loop
	//alpha = alpha-1;
	//alpha = alpha/50;
	//alpha = alpha-1;
	//alpha = pow(10,alpha);
	

	//dynamically allocate depending on cellsize
	bool charray[8];
	for(int i=0; i<cellsize; i++) {
		for(int j=0; j<cellsize; j++) {
			for(int k=0; k<cellsize; k++) {
				charray[4*i + 2*j + k] = true; 
			}
		}
	}

	charray[1] = false;
	//charray[2] = false;
	//charray[3] = false;	
	//charray[4] = false;

	double realenergy = 0;
	double kenergy=0;	

	for(int i=0; i<cellsize; i++) {
		for(int j=0; j<cellsize; j++) {
			for(int k=0; k<cellsize; k++) {
				realenergy += realsum(i,j,k,alpha,recip_cut,real_cut,cellsize,charray);
				kenergy += recsum(i,j,k,alpha,recip_cut,real_cut,cellsize,charray);
			}
		}
	}

	//measure total charge
	// double totalch=0;
	// for(int q=0;q<cellsize*cellsize*cellsize;q++){
	// 	if(charray[q]) {totalch++;}
	// 	else {totalch--;}
	// }
	// cout << "total charge: " << totalch << '\n';
	
	//net charge term
	//double netchterm = -(M_PI/2)*totalch/alpha/alpha;

	double selfenergy = selfint(alpha,cellsize);

	//cout << "" << realenergy << ' ';
	//cout << "" << kenergy << ' ';
	//cout << "" << selfenergy << ' ';
	cout << "" << realenergy + kenergy + selfenergy << '\n';

	//cout << "Ewald: " << totalenergy + netchterm << '\n';
	return 0;
	}
