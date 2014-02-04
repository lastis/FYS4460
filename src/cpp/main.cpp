#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CPhys.h"

using namespace std;
using namespace CPhys;
using namespace physical::constants;

void	force();

int main(int argc, const char *argv[]) {
	// NOTES:
	// Nc cannot be less than 1
	long	seed 	= 123;
	double 	sigma 	= 3.405; 	//E-10
	double 	b 	= 5.260/sigma; 	//E-10
	int 	Nc 	= 8;
	double	mass 	= 39.948*amu;
	double 	T	= 100;
	double	e0	= 119.8*K_B;
	double	v0	= sqrt(e0/m);
	double 	ran = 0;
	double 	vel = 0;
	Matrix 	r = Lattice::getFCC(Nc,b);
	Matrix 	v = Matrix(r);
	double** pR = r.getArrayPointer();
	double** pV = v.getArrayPointer();
	double	atoms = 4*Nc*Nc*Nc;
	// Output
	int	width 		= 15;
	int 	precision 	= 10;


	// Initialize the velocities
	for (int i = 0; i < atoms; i++) {
		vel = sigma/v0*Random::gauss(seed);
		v(i,0) = vel;
		vel = sigma/v0*Random::gauss(seed);
		v(i,1) = vel;
		vel = sigma/v0*Random::gauss(seed);
		v(i,2) = vel;
	}

	double 	a[3];
	double 	rij1[3];
	double 	rij2[3];
	double	lHalf = boxsize/2;
	Vector	fVec = Vector(atoms);

	// The verlet integration
	for (int t = 0; t < T; t+dt) {
		//  Calculate the force on the atoms
		for (int i = 0; i < atoms; i++) {
			for (int j = i; j < atoms; j++) {
				force(a,ri,rj);
			}
		}
		double vhalf[3];
		for (int i = 0; i < atoms; i++) {
			vhalf[0] = pV[i][0] + 0.5*dt*a[i][0];
			vhalf[1] = pV[i][1] + 0.5*dt*a[i][1];
			vhalf[2] = pV[i][2] + 0.5*dt*a[i][2];

			pR[i][0] += vhalf[0]*dt;
			pR[i][1] += vhalf[1]*dt;
			pR[i][2] += vhalf[2]*dt;
			if(pR[i][0] > boxsize){
				     pR[i][0] -= boxsize;
			}
			else if(pR[i][0] < 0){
				pR[i][0] += boxsize;
			}
			if(pR[i][1] > boxsize){
				pR[i][1] -= boxsize;
			}
			else if(pR[i][1] < 0){
				pR[i][1] += boxsize;
			}
			if(pR[i][2] > boxsize){
				pR[i][2] -= boxsize;
			}
			else if(pR[i][2] < 0){
				pR[i][2] += boxsize;
			}
		}
		//Recalculate forces

	}

	ofstream outFile("fcc.xyz");
	outFile << atoms << endl;
	outFile << "Argon atoms in a fcc latice" << endl;
	for (int i = 0; i < atoms; i++) {
		outFile << setw(width) << "Ar";
		outFile << setw(width) << setprecision(precision) 
			<< r(i,0);
		outFile << setw(width) << setprecision(precision) 
			<< r(i,1);
		outFile << setw(width) << setprecision(precision) 
			<< r(i,2);
		// Velocity
		outFile << setw(width) << setprecision(precision) 
			<< v(i,0);
		outFile << setw(width) << setprecision(precision) 
			<< v(i,1);
		outFile << setw(width) << setprecision(precision) 
			<< v(i,2);
		outFile << endl;
	}
	outFile.close();
	return 0;
}

void force(double* a, double ri, double rj, double L){
	for (int i = 0; i < atoms; i++) {
		for (int j = i+1; j < atoms; j++) {
			rij1[0] = pR[j][0] - pR[i][0];
			rij1[1] = pR[j][1] - pR[i][1];
			rij1[2] = pR[j][2] - pR[i][2];

			double r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
			double r2i = 1/r2;
			double r6i = r2i*r2i*r2i;
			double r12i = r6i*r6i;
		}
	}
};
