#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CPhys.h"

using namespace std;
using namespace CPhys;
using namespace physical::constants;

void force      (double**       a, double**     r, double     atoms, 
	         double   L);

void verlet     (double**       r, double**     v, double**       a, 
	         double        dt, int      atoms, double   L);

void printToFile(double**       r, double**     v, int        atoms, 
		 int 	    frame, int      width, int    precision);

Cube box;
double*	pointer;
double 	boxSize;
int 	boxes;
	
int main(int argc, const char *argv[]) {
	// NOTES:
	// Nc cannot be less than 1
	long	seed 	= 123;
	// Runtime
	double	dt	= 0.02;
	double	finalT  = 10*dt;
	// Simulation variables
	double 	sigma 	= 3.405; 	//E-10
	double 	b 	= 5.260/sigma; 	//E-10
	double	mass 	= 39.948*amu;
	double 	T	= 100;
	double	e0	= 119.8*K_B;
	double	v0	= sqrt(e0/m);
	double	stdDev 	= sqrt(K_B*T/m);
	// Atoms
	int 	Nc 	= 14;
	int	atoms 	= 4*Nc*Nc*Nc;
	double	L = b*Nc;
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
	// Allocate memory for box array
	boxSize = 3*sigma;
	boxes = L/boxSize + 1;
	box = Cube(boxes,boxes,boxes);

	pointer = Vector(atoms).getArrayPointer();

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
	force(pA,pR,atoms,L);
	// The verlet integration
	for (double i = 1; i*dt < finalT; i++) {
		verlet(pR,pV,pA,dt,atoms,L);
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

void refreshBoxes(double** r, int atoms, double L){
	// Reset the box, we want default value to be -1
	// (Note: we do not need to reset the pointer vector
	// because how the add method is created)
	box = -1;
	// Add
	int x = 0; int y = 0; int z = 0;
	for (int i = 0; i < atoms; i++) {
		// Division to find the box, then do cutoff to int
		x = r[i][0]/boxSize;
		y = r[i][1]/boxSize;
		z = r[i][2]/boxSize;
		pointer[i] = box(x,y,z);
		box(x,y,z) = i;
	}
};

void verlet(double**  r, double**     v, double**       a, 
	    double   dt, int      atoms, double   L){
	double vHalf[3] = {0,0,0};
	// Loop over all particles and update a vHalf speed.
	for (int i = 0; i < atoms; i++) {
		vHalf[0] = v[i][0] + 0.5*dt*a[i][0];
		vHalf[1] = v[i][1] + 0.5*dt*a[i][1];
		vHalf[2] = v[i][2] + 0.5*dt*a[i][2];

		r[i][0] += vHalf[0]*dt;
		r[i][1] += vHalf[1]*dt;
		r[i][2] += vHalf[2]*dt;

		// Periodic boundry conditions
		if(r[i][0] > L){
			r[i][0] -= L;
		}
		else if(r[i][0] < 0){
			r[i][0] += L;
		}
		if(r[i][1] > L){
			r[i][1] -= L;
		}
		else if(r[i][1] < 0){
			r[i][1] += L;
		}
		if(r[i][2] > L){
			r[i][2] -= L;
		}
		else if(r[i][2] < 0){
			r[i][2] += L;
		}
	}

	// Recalculate forces
	force(a,r,atoms,L);
	// Calculate the new velocities
	for (int i = 0; i < atoms; i++) {
		v[i][0] = vHalf[0] + 0.5*dt*a[i][0];
		v[i][1] = vHalf[1] + 0.5*dt*a[i][1];
		v[i][2] = vHalf[2] + 0.5*dt*a[i][2];
	}
};


void force1(double** a, double** r, int i, int j, double L){
	double rij[3] = {0,0,0};
	double lHalf = L/2.0;
	// Create the relative vector
	rij[0] = r[j][0] - r[i][0];
	rij[1] = r[j][1] - r[i][1];
	rij[2] = r[j][2] - r[i][2];

	// Check for closest path in periodic boundries
	if(rij[0] > lHalf){
		rij[0] -= L;
	}
	else if(rij[0] < -lHalf){
		rij[0] += L;
	}
	if(rij[1] > lHalf){
		rij[1] -= L;
	}
	else if(rij[1] < -lHalf){
		rij[1] += L;
	}
	if(rij[2] > lHalf){
		rij[2] -= L;
	}
	else if(rij[2] < -lHalf){
		rij[2] += L;
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
};

void force(double** a, double** r, double atoms, double L){
	refreshBoxes(r,atoms,L);
	int atom1 = 0; int atom2 = 0;
	int x = 0; int y = 0; int z = 0; 
	for (int i = 0; i < boxes; i++) {
	  for (int j = 0; j < boxes; j++) {
	    for (int k = 0; k < boxes; k++) {
	      // We choose box(i,j,k) and itterate through
	      // the boxes around this box
	      for (int x = i-1; x <= i+1; x++) {
	      	for (int y = j-1; y <= j+1; y++) {
	      	  for (int z = k-1; z <= k+1; z++) {
			// This is the first index of the atom
			// in the box(x,y,z)
			atom1 = box(i,j,k);
			atom2 = box(x,y,z);
			while (atom1 != -1) {
				while (atom2 != -1) {
					// Skip to next when the particle find
					// itself in the box
					if (atom1 == atom2){
						atom2 = pointer[atom2];
						continue;
					}	
					force1(a,r,atom1,atom2,L);
					// Done, now change j to the next 
					// particle in the box
					atom2 = pointer[atom2];
				}
				// Done, now change j to the next 
				// particle in the box
				atom1 = pointer[atom1];
			}
	      	  }
		}
	      }
	    }
	  }
	}

	/*

	int j = 0;
	int xCur = 0; int yCur = 0; int zCur = 0; 
	for (int i = 0; i < atoms; i++) {
		// Find the box which contains particle i
		xCur = r[i][0]/boxSize;
		yCur = r[i][1]/boxSize;
		zCur = r[i][2]/boxSize;
		// Itterate through the boxes on all sides
		for (int x = xCur-1; x <= xCur+1; x++) {
			for (int y = yCur-1; y <= yCur+1; y++) {
				for (int z = zCur-1; z <= zCur; z++) {
					if(z == -1) z = 0;
					// This is the first index of the atom
					// in the box(x,y,z)
					j = box[x][y][z];
					while (j != -1) {
						// Skip to next when the particle find
						// itself in the box
						if (j == i){
							j = pointer[j];
							continue;
						}	

						force1(a,r,i,j,L);

						// Done, now change j to the next 
						// particle in the box
						j = pointer[j];
					}
					// z loop ends here
				}
				// y loop ends here
			}
			// x loop ends here
		}
	}
	*/
};
