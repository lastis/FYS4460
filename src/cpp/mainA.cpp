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
	int Nx = 3;
	int Ny = 3;
	int Nz = 3;
	int cubes = (Nx-1)*(Ny-1)*(Nz-1);
	int bAtomsYZ = (Nz-1)*Nx + (Nz-1)*(Nx-1);
	int bAtomsXZ = (Ny-1)*Nz + (Ny-1)*(Nz-1);
	int bAtomsXY = (Nx-1)*Ny + (Nx-1)*(Ny-1);

	int atoms = cubes*4 + bAtomsXY + bAtomsXZ + bAtomsYZ + 1;
	Matrix r = Matrix(atoms,3);
	// Initiate r vector
	int cubeIndex = 0;
	for (int k = 0; k < Nz-1; k++) {
		for (int j = 0; j < Ny-1; j++) {
			for (int i = 0; i < Nx-1; i++) {
				cubeIndex = k*((Ny-1)*(Nx-1))+j*(Nx-1)+i;
				r(cubeIndex,0) = i*b;
				r(cubeIndex,1) = j*b;
				r(cubeIndex,2) = k*b;

				r(cubeIndex+1,0) = i*b + 0.5*b;
				r(cubeIndex+1,1) = j*b + 0.5*b;
				r(cubeIndex+1,2) = k*b;

				r(cubeIndex+2,0) = i*b;
				r(cubeIndex+2,1) = j*b + 0.5*b;
				r(cubeIndex+2,2) = k*b + 0.5*b;

				r(cubeIndex+3,0) = i*b + 0.5*b;
				r(cubeIndex+3,1) = j*b;
				r(cubeIndex+3,2) = k*b + 0.5*b;
			}
		}
	}
	int index = cubes*4;
	for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			// TODO continue, this is where i stopped
			r(index+)
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
