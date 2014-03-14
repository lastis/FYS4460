using namespace physical::constants;
using namespace physical::unit;

/////////////////////////////////////////////////
// 		State Variables
/////////////////////////////////////////////////

const double 	sigma 	= 3.405; 	//E-10
const double 	b 	= 5.720/sigma; 	//E-10
const double 	mass 	= 39.948*amu;
const double 	e0	= 119.8*K_B;
const double	T0	= e0/K_B;
const double 	T	= 0.851*T0;
const double 	v0	= sqrt(e0/mass);
const double 	stdDev 	= sqrt(K_B*T/mass);
const int	Nc	= 20;

/////////////////////////////////////////////////
// 		Simulation Variables
/////////////////////////////////////////////////

const double 	dt	= 0.02;
const double 	finalT  = 1000*dt;
const int 	dumpRate 	= 100;
//const double	targetT = 0.851*T0;
const double	targetT = 1.20*T0;
const double	tau	= 20*dt;
const double    boxSize = 2;

const bool	usePores = false;
const double	poreRadius = 1*nm/(sigma*1e-10);
const int	poreCnt = 1;

const bool	useCylinders = true;
const double	poreCylRadius = 2*nm/(sigma*1e-10);
const int	poreCylCnt = 4;
