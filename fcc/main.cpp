#include "ewald.h"

using namespace std;

int main(int argc, char *argv[]){

	time_t start = time(NULL);

	double alpha;
	int real_cut, recip_cut, cellsize;
	
	cellsize = atoi(argv[1]);
	alpha = strtod(argv[2],NULL);
	real_cut = atoi(argv[3]);
	recip_cut = atoi(argv[4]);

	int bsize = 4;
	int fsize = 4;

	int N = fsize*bsize*cellsize*cellsize*cellsize;
	
	double off[12] = {0, 0, 0,
			0, 0.25, 0.25,
			0.25, 0, 0.25,
			0.25, 0.25, 0}; //not necessary here?

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
			for(int p=0; p<fsize; p++){
				realenergy += realsum(i,j,k,m,p,alpha,real_cut,cellsize,bsize,&totNNenergy,intmat);
				kenergy += recsum(i,j,k,m,p,alpha,recip_cut,cellsize,bsize,intmat);
			}
			}
			}
		}
	}

	double selfenergy = selfint(alpha,cellsize,bsize);

	cout << "NN: " << totNNenergy << '\n';
	cout << "selfenergy: " << selfenergy << '\n';
	//cout << "Ewald w/o selfen: " << realenergy + kenergy << '\n';
	cout << "Ewald : " << realenergy + kenergy + selfenergy << '\n';


	//write to a binary file
	FILE* bmatstream;
	char bmatname[50];
	sprintf(bmatname,"IntmatJ0_a%d_r%d_k%d.bin",cellsize,real_cut,recip_cut);
	
	bmatstream = fopen(bmatname, "wb");
	fwrite(intmat, sizeof(double), N*N, bmatstream);
	fclose(bmatstream);

	/*char matname[50];
	sprintf(matname,"Intmat_a%d_r%d_k%d.txt",cellsize,real_cut,recip_cut);

	ofstream matfile;
	matfile.open(matname);
	for(int i=0;i<N*N;i++) {
		matfile << intmat[i] << '\n';
	}
	matfile.close();*/

	delete[] intmat;

	time_t end = time(NULL);

	printf("cellsize: %d\n",cellsize);
	printf("alpha: %.2f\n",alpha);

	printf("\n");
	printf("time: %5.3f\n", (double)(end - start));
	printf("in hours: %2.2f\n", (double)(end - start)/3600);

	return 0;
}
