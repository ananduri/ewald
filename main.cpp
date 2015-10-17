#include "ewald.h"

using namespace std;

int main(int argc, char *argv[]){

	double alpha;
	int real_cut, recip_cut, cellsize;
	
	cellsize = atoi(argv[1]);
	alpha = strtod(argv[2],NULL);
	real_cut = atoi(argv[3]);
	recip_cut = atoi(argv[4]);

	int bsize = 4;

	int N = bsize*cellsize*cellsize*cellsize;
	
	bool* charray = NULL;
	charray = new bool[N];

	for(int i=0;i<N;i++) {
		charray[i] = true;
	}
	


	//for(int i =0;i<N;i++){
	//	cout << charray[i] << endl;
	//}


	//charray[0] = true;
	//charray[20] = true;

	/*for(int i=0; i<cellsize; i++) {
		for(int j=0; j<cellsize; j++) {
			for(int k=0; k<cellsize; k++) {
				for (int l=0; l<bsize; l++) {
				cout << (charray[bsize*cellsize*cellsize*i + bsize*cellsize*j + bsize*k + l] ? "true" : "false") << endl;
			}
			}
		}
	}*/

	
	double* intmat = NULL;
	intmat = new double[N*N];
	for(int i=0;i<N*N;i++){
		intmat[i] = 0.0;
	}



	double realenergy = 0;
	double kenergy=0;	
	double totNNenergy=0;

	for(int i=0; i<cellsize; i++) {
		for(int j=0; j<cellsize; j++) {
			for(int k=0; k<cellsize; k++) {
			for(int m=0; m<bsize; m++){
				realenergy += realsum(i,j,k,m,alpha,real_cut,cellsize,bsize,charray,&totNNenergy,intmat);
				kenergy += recsum(i,j,k,m,alpha,recip_cut,cellsize,bsize,charray,intmat);
			}
			}
		}
	}
	delete[] charray;

	double selfenergy = selfint(alpha,cellsize,bsize);

	cout << "NN: " << totNNenergy << '\n';
	cout << "selfenergy: " << selfenergy << '\n';
	cout << "Ewald w/o selfen: " << realenergy + kenergy << '\n';
	cout << "Ewald : " << realenergy + kenergy + selfenergy << '\n';


	//write to a binary file in first piece
	FILE* bmatstream;
	char bmatname[50];
	sprintf(bmatname,"Intmat_a%d_r%d_k%d.bin",cellsize,real_cut,recip_cut);
	
	bmatstream = fopen(bmatname, "wb");
	fwrite(intmat, sizeof(double), N*N, bmatstream);
	fclose(bmatstream);



	char matname[50];
	sprintf(matname,"Intmat_a%d_r%d_k%d.txt",cellsize,real_cut,recip_cut);

	ofstream matfile;
	matfile.open(matname);
	for(int i=0;i<N*N;i++) {
		matfile << intmat[i] << '\n';
	}
	matfile.close();

	delete[] intmat;
	return 0;
}
