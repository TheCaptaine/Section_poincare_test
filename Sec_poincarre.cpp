/* COMPILATION:

gfortran -Wall -c radau.f90
g++ -c Sec_poincarre.cpp
g++ -Wall radau.o Sec_poincarre.o -o Sec_poincarre.x -lgfortran

*/

#include <iostream>
#include <math.h>

extern "C" {
   void RA15_I(int, int, double[], double*, void(int,double[],double*,double[]), double, double, void(int,double[],double*));
}

double alea_uniform(double, double);
void force(int, double[], double*, double[]);
void force2(int, double[], double*, double[]);
double hamiltonian(double[]);
double hamiltonian2(double[]);
void printFunc(int, double[], double*);
void silentprintFunc(int, double[], double*);

using namespace std;

int Nstep;
int Nsave = 1;

int main(int argc, char *argv[]){
	srand(1234);
	const int dimX = 4;
	double X0[dimX], t0, H0;
	double X[dimX], t;

	int LL;
	double dTini, Tspan, E;

	int k;

	printf("# t q1 q2 p1 p2\n");

	dTini = 0.1;
	Tspan = 1.e3;

	LL = 10;

	t0 = 0.0;
	E = strtod(argv[1], NULL);
	for (int i = 0; i < 1000; i++) {
		printf("# %d\n", i+1);
		do {
			X0[0] = alea_uniform(0.15, 0.4);
			X0[1] = alea_uniform(0.15, 0.4);
			X0[2] = alea_uniform(0.15, 0.4);
   		X0[3] = 2*(E-pow(X0[0], 2)*X0[1]+1./3*pow(X0[1], 3))-(pow(X0[0], 2)+pow(X0[1], 2)+pow(X0[2], 2));
   	} while(X0[3] < 0);
   	X0[3] = pow(X0[3], 0.5);
   /*printf("# %22.15e\n", X0[3]);
   X0[3] = pow(X0[3], 0.5);
   printf("# %22.15e\n", X0[3]);
   H0 = hamiltonian2(X0);
	printf("# %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", E, H0, X0[0], X0[1], X0[2], X0[3]);*/

		t = t0;
		X[0] = X0[0];
		X[1] = X0[1];
		X[2] = X0[2];
		X[3] = X0[3];
		Nstep = 0;


		printFunc(dimX, X, &t);
		RA15_I(LL, dimX, X, &t, force2, dTini, Tspan, printFunc);
	}
	return 0;
}

double alea_uniform(double a, double b) {
	double x = rand()/(RAND_MAX+1.0);
	return (b-a)*x+a;
}

void force(int dimX, double X[], double* t, double F[]){
   	F[0] = X[1];
   	F[1] = -sin(X[0]);
}

void force2(int dimX, double X[], double* t, double F[]){
   	F[0] = X[2];
   	F[1] = X[3];
   	F[2] = -(X[0]+2*X[0]*X[1]);
   	F[3] = -(X[1]+pow(X[0], 2)-pow(X[1], 2));
}

double hamiltonian(double X[]){
   	return 0.5*X[1]*X[1] - cos(X[0]);
}

double hamiltonian2(double X[]){
   	return 0.5*(pow(X[0], 2)+pow(X[1], 2)+pow(X[2], 2)+pow(X[3], 2))+pow(X[0], 2)*X[1]-1./3*pow(X[1], 3);
}

void printFunc(int dimX, double X[], double* t){
   	if(Nstep%Nsave == 0){
      		printf("%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", *t, X[0], X[1], X[2], X[3], hamiltonian2(X));
   	}
   	Nstep++;
}

void silentprintFunc(int dimX, double X[], double* t){
   	Nstep++;
}
