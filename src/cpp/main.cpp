#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CPhys.h"

using namespace std;
using namespace CPhys;
using namespace physical::constants;

void force      (double**       a, double**     r, double     atoms, 
	         double   boxSize);

void verlet     (double**       r, double**     v, double**       a, 
	         double        dt, int      atoms, double   boxSize);

void printToFile(double**       r, double**     v, int        atoms, 
		 int 	    frame, int      width, int    precision);

int main(int argc, const char *argv[]) {
	// NOTES:
	// Nc cannot be less than 1
	long	seed 	= 123;
	// Runtime
	double	dt	= 0.02;
	double	finalT  = 100*dt;
	// Simulation variables
	double 	sigma 	= 3.405; 	//E-10
	double 	b 	= 5.260/sigma; 	//E-10
	double	mass 	= 39.948*amu;
	double 	T	= 100;
	double	e0	= 119.8*K_B;
	double	v0	= sqrt(e0/m);
	double	stdDev 	= sqrt(K_B*T/m);
	// Atoms
	int 	Nc 	= 8;
	int	atoms 	= 4*Nc*Nc*Nc;
	double	boxSize = b*Nc;
	// Vectors
	Matrix 	r = Lattice::getFCC(Nc,b); // r = n x m : n = atoms : m = 3
	Matrix 	v = Matrix(atoms,3);
	Matrix 	a = Matrix(atoms,3);
	double** pR = r.getArrayPointer();
	double** pV = v.getArrayPointer();
	double** pA = a.getArrayPointer();
	// Output
	int	width 		= 16;
	int 	precision 	= 8;

	// Initialize the velocities
	double 	vel = 0;
	for (int i = 0; i < atoms; i++) {
		vel = stdDev/v0*Random::gauss(seed);
		v(i,0) = vel;
		vel = stdDev/v0*Random::gauss(seed);
		v(i,1) = vel;
		vel = stdDev/v0*Random::gauss(seed);
		v(i,2) = vel;
	}

	// Do first printout of first state
	printToFile(pR,pV,atoms,0,width,precision);
	// Calculate the total force on all atoms once first
	force(pA,pR,atoms,boxSize);
	// The verlet integration
	for (double i = 1; i*dt < finalT; i++) {
		verlet(pR,pV,pA,dt,atoms,boxSize);
		printToFile(pR,pV,atoms,i,width,precision);
	}

	return 0;
}

void printToFile(double**     r, double**     v, int     atoms, 
		 int 	  frame, int      width, int precision){
	char filename[64];
	sprintf(filename,"./../../res/result/%04d.xyz",frame);

	ofstream outFile(filename);
	outFile << atoms << endl;
	outFile << "Argon atoms in a fcc latice" << endl;
	for (int i = 0; i < atoms; i++) {
		outFile << "Ar";
		outFile << setw(width) << setprecision(precision) 
			<< r[i][0];
		outFile << setw(width) << setprecision(precision) 
			<< r[i][1];
		outFile << setw(width) << setprecision(precision) 
			<< r[i][2];
		// Velocity
		outFile << setw(width) << setprecision(precision) 
			<< v[i][0];
		outFile << setw(width) << setprecision(precision) 
			<< v[i][1];
		outFile << setw(width) << setprecision(precision) 
			<< v[i][2];
		outFile << endl;
	}
	outFile.close();
};

void verlet(double**  r, double**     v, double**       a, 
	     double   dt, int      atoms, double   boxSize){
	// Loop over all particles and update a vHalf speed.
	double vHalf[3] = {0,0,0};
		for (int i = 0; i < atoms; i++) {
			vHalf[0] = v[i][0] + 0.5*dt*a[i][0];
			vHalf[1] = v[i][1] + 0.5*dt*a[i][1];
			vHalf[2] = v[i][2] + 0.5*dt*a[i][2];

			r[i][0] += vHalf[0]*dt;
			r[i][1] += vHalf[1]*dt;
			r[i][2] += vHalf[2]*dt;

			// Periodic boundry conditions
			if(r[i][0] > boxSize){
				     r[i][0] -= boxSize;
			}
			else if(r[i][0] < 0){
				r[i][0] += boxSize;
			}
			if(r[i][1] > boxSize){
				r[i][1] -= boxSize;
			}
			else if(r[i][1] < 0){
				r[i][1] += boxSize;
			}
			if(r[i][2] > boxSize){
				r[i][2] -= boxSize;
			}
			else if(r[i][2] < 0){
				r[i][2] += boxSize;
			}
		}

		// Recalculate forces
		force(a,r,atoms,boxSize);
		// Calculate the new velocities
		for (int i = 0; i < atoms; i++) {
			v[i][0] = vHalf[0] + 0.5*dt*a[i][0];
			v[i][1] = vHalf[1] + 0.5*dt*a[i][1];
			v[i][2] = vHalf[2] + 0.5*dt*a[i][2];
		}
};

void force(double** a, double** r, double atoms, double boxSize){
	double rij[3] = {0,0,0};
	double lHalf = boxSize/2.0;
	for (int i = 0; i < atoms; i++) {
		for (int j = i+1; j < atoms; j++) {
			// Create the relative vector
			rij[0] = r[j][0] - r[i][0];
			rij[1] = r[j][1] - r[i][1];
			rij[2] = r[j][2] - r[i][2];

			// Check for closest path in periodic boundries
			if(rij[0] > lHalf){
				rij[0] -= boxSize;
			}
		    	else if(rij[0] < -lHalf){
				rij[0] += boxSize;
			}
			if(rij[1] > lHalf){
				rij[1] -= boxSize;
			}
		    	else if(rij[1] < -lHalf){
				rij[1] += boxSize;
			}
			if(rij[2] > lHalf){
				rij[2] -= boxSize;
		    	}
		    	else if(rij[2] < -lHalf){
			    	rij[2] += boxSize;
			}

			double r2 = rij[0]*rij[0] 
				  + rij[1]*rij[1] 
				  + rij[2]*rij[2];
			double r2i = 1/r2;
			double r6i = r2i*r2i*r2i;
			double r12i = r6i*r6i;
			a[i][0] -= 24*(2*r12i-r6i) * r2i * rij[0];
			a[i][1] -= 24*(2*r12i-r6i) * r2i * rij[1];
			a[i][2] -= 24*(2*r12i-r6i) * r2i * rij[2];

			a[j][0] += 24*(2*r12i-r6i) * r2i * rij[0];
			a[j][1] += 24*(2*r12i-r6i) * r2i * rij[1];
			a[j][2] += 24*(2*r12i-r6i) * r2i * rij[2];
		}
	}
};
