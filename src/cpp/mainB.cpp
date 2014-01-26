#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CPhys.h"

using namespace std;
using namespace CPhys;
using namespace physical::constants;

int main(int argc, const char *argv[])
{
	// NOTES:
	// Nc cannot be less than 1
	// Counting is done by adding all atoms from all unit cells
	// then adding the sides except the edges, then add the 
	// edges (they overlap on the last atom)
	long	seed 	= 123;
	double 	b 	= 5.260;
	int 	Nc 	= 3;
	int 	cubes 	= (Nc-1)*(Nc-1)*(Nc-1);
	int 	atomsSide = 2*3*(Nc-1)*(Nc-1);
	int 	atoms 	= cubes*4 + atomsSide + 3*Nc - 2;
	double	mass 	= 39.948*amu;
	double 	sigma 	= sqrt(K_B*100/mass);
	double 	ran = 0;
	double 	vel = 0;
	// Output
	int	width 		= 15;
	int 	precision 	= 10;

	Matrix 	r = Matrix(atoms,3);
	Matrix 	v = Matrix(atoms,3);
	// Initiate r vector
	int index = 0;
	for (int k = 0; k < Nc-1; k++) {
		for (int j = 0; j < Nc-1; j++) {
			for (int i = 0; i < Nc-1; i++) {
				r(index,0) = i*b;
				r(index,1) = j*b;
				r(index,2) = k*b;
				index++;

				r(index,0) = i*b + 0.5*b;
				r(index,1) = j*b + 0.5*b;
				r(index,2) = k*b;
				index++;

				r(index,0) = i*b;
				r(index,1) = j*b + 0.5*b;
				r(index,2) = k*b + 0.5*b;
				index++;

				r(index,0) = i*b + 0.5*b;
				r(index,1) = j*b;
				r(index,2) = k*b + 0.5*b;
				index++;
			}
		}
	}
	
	// Add position of atoms on the sides
	for (int x = 0; x < Nc-1; x++) {
		for (int y = 0; y < Nc-1; y++) {
			r(index,0) = (Nc-1)*b;
			r(index,1) =  x*b;
			r(index,2) =  y*b;
			index++;
			
			r(index,0) = (Nc-1)*b;
			r(index,1) =  x*b + 0.5*b;
			r(index,2) =  y*b + 0.5*b;
			index++;

			r(index,0) =  x*b;
			r(index,1) = (Nc-1)*b;
			r(index,2) =  y*b;
			index++;
			
			r(index,0) =  x*b + 0.5*b;
			r(index,1) = (Nc-1)*b;
			r(index,2) =  y*b + 0.5*b;
			index++;

			r(index,0) =  x*b;
			r(index,1) =  y*b;
			r(index,2) = (Nc-1)*b;
			index++;
			
			r(index,0) =  x*b + 0.5*b;
			r(index,1) =  y*b + 0.5*b;
			r(index,2) = (Nc-1)*b;
			index++;
		}
	}

	// Add position of atoms on the edges
	for (int a = 0; a < Nc-1; a++) {
		r(index, 0) = (Nc-1)*b;
		r(index, 1) = a*b;
		r(index, 2) = (Nc-1)*b;
		index++;
		
		r(index, 0) = a*b;
		r(index, 1) = (Nc-1)*b;
		r(index, 2) = (Nc-1)*b;
		index++;
		
		r(index, 0) = (Nc-1)*b;
		r(index, 1) = (Nc-1)*b;
		r(index, 2) = a*b;
		index++;
	}
	// Add the last atom on the crossing of the edges
	r(index,0) = (Nc-1)*b;
	r(index,1) = (Nc-1)*b;
	r(index,2) = (Nc-1)*b;

	for (int i = 0; i < atoms; i++) {
		vel = sigma*Random::gauss(seed);
		v(i,0) = vel;
		vel = sigma*Random::gauss(seed);
		v(i,1) = vel;
		vel = sigma*Random::gauss(seed);
		v(i,2) = vel;
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
