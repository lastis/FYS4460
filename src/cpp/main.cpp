#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CPhys.h"
#include "Parameters.h"

using namespace std;
using namespace CPhys;

void force      (double**       a, double**     r, Cube& box, double*& list);

void verlet     (double**       r, double**     v, double**       a, 
	         Cube& box, double*& list);

void printToFile(double**       r, double**     v,
		 int 	    frame, int      width, int    precision);
void initState(Matrix& r, Matrix& v, Matrix& a, long seed);
void initBoxes(Cube& box, double*& list);
void applyForce(double** a, double** r, int i, int j);
void refreshBoxes(double** r, Cube& box, double*& list);
///////////////////////////////////////////
// 		Read only variables
///////////////////////////////////////////

// State variables
extern const double 	sigma;
extern const double 	b;
extern const double	mass;
extern const double 	T;
extern const double	e0;
extern const double	v0;
extern const double	stdDev;
// Simulation variables
extern const double 	dt;
extern const double 	finalT;
extern const int 	Nc;
extern const int 	atoms;
extern const double 	L;
extern const double 	boxSize;
extern const int 	boxes;

int main(int argc, const char *argv[]) {
	long	seed 	= 123;
	// Vectors
	Matrix 	r; double** ppR;
	Matrix 	v; double** ppV;
	Matrix 	a; double** ppA;
	Cube 	box;
	double* list;
	// Output
	int	width 		= 16;
	int 	precision 	= 8;

	// Initialize state
	initState(r,v,a,seed);
	initBoxes(box,list);

	// Get the underlying list for speed. 
	ppR = r.getArrayPointer();
	ppV = v.getArrayPointer();
	ppA = a.getArrayPointer();

	///////////////////////////////////////////
	// 		Simulation
	///////////////////////////////////////////
	
	printToFile(ppR,ppV,0,width,precision);
	// Calculate the total force on all atoms once first
	force(ppA,ppR,box,list);
	// Verlet integration
	for (double i = 1; i*dt < finalT; i++) {
		verlet(ppR,ppV,ppA,box,list);
		printToFile(ppR,ppV,i,width,precision);
	}
	return 0;
}

void initState(Matrix& r, Matrix& v, Matrix& a, long seed){
	r = Lattice::getFCC(Nc,b); // r = n x m : n = atoms : m = 3
	v = Matrix(atoms,3);
	a = Matrix(atoms,3);

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
	// Remove linear moment
	double **ppV = v.getArrayPointer();
	double usum = 0; double vsum = 0; double wsum = 0;
	for(int i = 0; i < atoms; i++){
		usum += ppV[i][0];
		vsum += ppV[i][1];
		wsum += ppV[i][2];
	}  
	double du = usum/atoms; double dv = vsum/atoms; double dw = wsum/atoms;
	for(int i = 0; i < atoms; i++){
		ppV[i][0] -= du;
		ppV[i][1] -= dv;
		ppV[i][2] -= dw;
	}
}

void initBoxes(Cube& box, double*& list){
	box 	= Cube(boxes,boxes,boxes);
	//list 	= Vector(atoms).getArrayPointer();
	list 	= new double[atoms];
}

void printToFile(double**     r, double**     v,
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

void refreshBoxes(double** r, Cube& box, double*& list){
	// Reset the box, we want default value to be -1
	// (Note: we do not need to reset the list vector
	// because how the add method is created)
	box = -1;
	// Add
	int x = 0; int y = 0; int z = 0;
	for (int i = 0; i < atoms; i++) {
		// Division to find the box, then do cutoff to int
		x = r[i][0]/boxSize;
		y = r[i][1]/boxSize;
		z = r[i][2]/boxSize;
		list[i] = box(x,y,z);
		box(x,y,z) = i;
	}
}

void verlet(double**  r, double**     v, double**       a, Cube& box, double*& list){
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
	force(a,r,box,list);
	// Calculate the new velocities
	for (int i = 0; i < atoms; i++) {
		v[i][0] = vHalf[0] + 0.5*dt*a[i][0];
		v[i][1] = vHalf[1] + 0.5*dt*a[i][1];
		v[i][2] = vHalf[2] + 0.5*dt*a[i][2];
	}
}


void applyForce(double** a, double** r, int i, int j){
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

}

void force(double** a, double** r, Cube& box, double*& list){
	refreshBoxes(r,box,list);
	int atom1 = 0; int atom2 = 0;
	int x = 0; int y = 0; int z = 0; 
	int i = 0; int j = 0; int k = 0;

	for (int atom1 = 0; atom1 < atoms; atom1++) {
		// Find the box of "atom1"
		i = r[atom1][0];
		j = r[atom1][1];
		k = r[atom1][2];

		for (int x = i-1; x <= i+1; x++) {
			for (int y = j-1; y <= j+1; y++) {
				for (int z = k-1; z <= k+1; z++) {
					// This is the first index of the atom2
					// in the box(x,y,z)
					atom2 = box(x,y,z);
					while (atom2 != -1) {
						// Skip to next when the particle find
						// itself in the box
						if (atom1 == atom2){
							atom2 = list[atom2];
							continue;
						}	
						applyForce(a,r,atom1,atom2);

						// Done, now change j to the next 
						// particle in the box
						atom2 = list[atom2];
					} // Second while loop ends here

				} // z loop ends here
			} // y loop ends here
		} // x loop ends here
	}
	/*
	refreshBoxes(r,box,list);
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
			  			atom2 = list[atom2];
			  			continue;
		        		}	
		        		applyForce(a,r,atom1,atom2);

		        		// Done, now change j to the next 
		        		// particle in the box
		        		atom2 = list[atom2];
		      		} // Second while loop ends here

		      		// Done, now change j to the next 
		      		// particle in the box
		      		atom1 = list[atom1];

		    	} // First while loop ends here

	      	  } // z loop ends here
		} // y loop ends here
	      } // x loop ends here
	    } // k loop ends here
	  } // j loop ends here
	} // i loop ends here
	*/
}
