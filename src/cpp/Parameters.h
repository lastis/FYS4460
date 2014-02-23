//using namespace physical::constants;
// State variables
const double    amu = 1.6605402e-27; 
const double    K_B = 1.380658e-23; 
const double 	sigma 	= 3.405; 	//E-10
const double 	b 	= 10;//10000*5.260/sigma; 	//E-10
const double 	mass 	= 39.948*amu;
const double 	T	= 100;
const double 	e0	= 119.8*K_B;
const double 	v0	= sqrt(e0/mass);
const double 	stdDev 	= sqrt(K_B*T/mass);
// Simulation variables
const double 	dt	= 0.02;
const double 	finalT  = 100*dt;
int     	    Nc 	= 1;
int      	    natoms 	= 4*Nc*Nc*Nc;
double       	L 	= b*Nc;
//const double 	boxSize = 3*sigma;
//const int 	boxes 	= L/boxSize + 1;
double 	        boxSize = L/2;
int 	        boxes 	= L/boxSize + 1;
