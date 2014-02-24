using namespace physical::constants;
// State variables
const double 	sigma 	= 3.405; 	//E-10
const double 	b 	= 10;//5.260/sigma; 	//E-10
const double 	mass 	= 39.948*amu;
const double 	T	= 100;
const double 	e0	= 119.8*K_B;
const double 	v0	= sqrt(e0/mass);
const double 	stdDev 	= sqrt(K_B*T/mass);
// Simulation variables
const double 	dt	= 0.02;
const double 	finalT  = 100*dt;
const int 	dumpRate 	= 10;

