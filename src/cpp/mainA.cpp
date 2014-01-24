#include <fstream>
#include <iostream>
#include <iomanip>
#include "CPhys.h"

using namespace std;
using namespace CPhys;

int main(int argc, const char *argv[])
{
	double b = 1;
	// Nx, Ny and Nz cannot be 1 or less
	int Nc = 2;
	int cubes = (Nc-1)*(Nc-1)*(Nc-1);
	// Three sides with pairs of atoms
	int atomsSide = 2*3*(Nc-1)*(Nc-1);

	// Atoms from all the unit cells, then add the sides except 
	// the edges, then add the edges (they overlap on the last atom)
	int atoms = cubes*4 + atomsSide + 3*Nc - 2;

	Matrix r = Matrix(atoms,3);
	// Initiate r vector
	int index = 0;
	for (int k = 0; k < Nc-1; k++) {
		for (int j = 0; j < Nc-1; j++) {
			for (int i = 0; i < Nc-1; i++) {
				//cubeIndex = k*((Nc-1)*(Nc-1))+j*(Nc-1)+i;
				//cubeIndex = 4*cubeIndex;

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
				index++
			}
		}
	}
	// TODO I stopped here, I changed the limits to Nc
	// and I need to merge these three loops. 
	
	// Add position of atoms on the sides
	for (int j = 0; j < Nc-1; j++) {
		for (int k = 0; k < Nc-1; k++) {
			r(index,0) = Nx*b;
			r(index,1) =  j*b;
			r(index,2) =  k*b;
			index++;
			
			r(index,0) = Nx*b;
			r(index,1) =  j*b + 0.5*b;
			r(index,2) =  k*b + 0.5*b;
			index++
		}
	}
	for (int i = 0; i < Nc-1; i++) {
		for (int k = 0; k < Nc-1; k++) {
			r(index,0) =  i*b;
			r(index,1) = Ny*b;
			r(index,2) =  k*b;
			index++
			
			r(index,0) =  i*b + 0.5*b;
			r(index,1) = Ny*b;
			r(index,2) =  k*b + 0.5*b;
			index++
		}
	}
	for (int i = 0; i < Nc-1; i++) {
		for (int j = 0; j < Nc-1; j++) {
			r(index,0) =  i*b;
			r(index,1) =  j*b;
			r(index,2) = Nz*b;
			
			r(index,0) =  i*b + 0.5*b;
			r(index,1) =  j*b + 0.5*b;
			r(index,2) = Nz*b;
		}
	}


	ofstream outFile("fcc.xyz");
	outFile << atoms << endl;
	outFile << "Argon atoms in a fcc latice" << endl;
	for (int i = 0; i < atoms; i++) {
		outFile << setw(5) << "Ar";
		outFile << setw(5) << setprecision(4) 
			<< r(i,0);
		outFile << setw(5) << setprecision(4) 
			<< r(i,1);
		outFile << setw(5) << setprecision(4) 
			<< r(i,2);
		// Velocity
		outFile << setw(5) << setprecision(4) 
			<< 0;
		outFile << setw(5) << setprecision(4) 
			<< 0;
		outFile << setw(5) << setprecision(4) 
			<< 0;
		outFile << endl;
	}
	outFile.close();
	return 0;
}
