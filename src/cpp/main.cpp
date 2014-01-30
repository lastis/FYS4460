#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CPhys.h"

using namespace std;
using namespace CPhys;
using namespace physical::constants;

int main(int argc, const char *argv[]) {
	// NOTES:
	// Nc cannot be less than 1
	long	seed 	= 123;
	double 	b 	= 5.260;
	int 	Nc 	= 8;
	double	mass 	= 39.948*amu;
	double 	sigma 	= sqrt(K_B*100/mass);
	double 	ran = 0;
	double 	vel = 0;
	double	e0	= 1;
	Matrix 	r = Lattice::getFCC(Nc,b);
	Matrix 	v = Matrix(r);
	double	atoms = r.getLength();
	// Output
	int	width 		= 15;
	int 	precision 	= 10;


	// Initialize the velocities
	for (int i = 0; i < atoms; i++) {
		vel = sigma*Random::gauss(seed);
		v(i,0) = vel;
		vel = sigma*Random::gauss(seed);
		v(i,1) = vel;
		vel = sigma*Random::gauss(seed);
		v(i,2) = vel;
	}

	/*
	// Verlet integration
	double ir = 0;
	double ir2 = 0;
	double ir6 = 0;
	double ir12 = 0;
	for (int i = 0; i < atoms; i++) {
		for (int j = 0; j < atoms; j++) {
			ir = 1/(r(j,0)+r(j,1)+r(j,2));
			ir2 = ir*ir;
			ir6 = ir2*ir2*ir2;
			ir12 = ir6*ir6;
			F = 24*(2*ir12-ir6);
			Fx = F*ir2*x;
			Fy = F*ir2*y;
			Fz = F*ir2*z;
		}
	}
	*/

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
