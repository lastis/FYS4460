using namespace physical::constants;
using namespace physical::unit;

/////////////////////////////////////////////////
// 		State Variables
/////////////////////////////////////////////////

const double 	sigma 	= 3.405; 	//E-10
const double 	b 	= 5.720/sigma; 	//E-10
const double 	mass 	= 39.948*amu;
const double 	T	= 1;
const double 	e0	= 119.8*K_B;
const double 	v0	= sqrt(e0/mass);
const double 	stdDev 	= sqrt(K_B*T/mass);
const int	Nc	= 8;

/////////////////////////////////////////////////
// 		Simulation Variables
/////////////////////////////////////////////////

const double 	dt	= 0.02;
const double 	finalT  = 200*dt;
const int 	dumpRate 	= 10;
const double	targetT = 0.851;
const double	tau	= 20*dt;
const double	poreRadius = 2*nm/(sigma*1e-10);
const int	poreCnt = 5;
