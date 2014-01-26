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
	long	seed 	= 123;
	double 	b 	= 5.260;
	int 	Nc 	= 8;
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

	Matrix 	r = Lattice::getFCC(Nc,b);
	Matrix 	v = Matrix(r);

	// Initialize the velocities
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
