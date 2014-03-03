#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include "CPhys.h"
#include "Parameters2.h"
using namespace std;
using namespace CPhys;

void doVerletIntegration();

void writeState(int frameNum);

void loadState(int frameNum);

void applyForce(int i, int j);

void computeForce();

void computeEnergy();

void computeTemperature(bool write);

void computeRadialDistr(double lHalf);

void computePressure();

void applyBerendsenThermostat();

void refreshBoxes();

void createPores();

void printSystemStats(int frameNum);

void writeFrozenState();

void initLocalParam();

void createCylindricPores();

/////////////////////////////////////////////////
// 		Simulation Variables
/////////////////////////////////////////////////

extern const double 	dt;
extern const double 	finalT;
extern const int 	dumpRate;
extern const double	targetT;
extern const double	tau;

extern const bool	usePores;
extern const double	poreRadius;
extern const int	poreCnt;

extern const bool	useCylinders;
extern const double	poreCylRadius;
extern const int	poreCylCnt;
/////////////////////////////////////////////////
// 		State Variables
/////////////////////////////////////////////////

extern const double 	sigma;
extern const double 	b;
extern const double 	mass;
extern const double 	T;
extern const double 	e0;
extern const double 	v0;
extern const double 	stdDev;
extern const int	Nc;


Cube	cBox;
Matrix	mState;
double** mpState;
double** mpA;
double  boxSize;
double 	L;
double	curT;
int 	natoms;
int 	boxes;
int*	vpLinkedList;
bool*	vpFrozenAtoms;

long	seed = 120;
double 	aij[3];
double 	rij[3];
double 	r2;
double 	r2i;
double 	r6i;
double 	r12i;
double 	sumPressure;
char 	tmpstr[64]; 

char outputPath[256];

int main(int nargs, char** argsv){
	/* Read starting point and initialize corresponding variables */
	int 	frameNum = atoi(argsv[1]);
    
    	// Set outputpath
    	sprintf(outputPath, "../../res/");
    	//sprintf(outputPath, "/home/andrenos/work/FYS4460/");
    
	loadState(frameNum); 
	initLocalParam();
	if (usePores) createPores();
	if (useCylinders) createCylindricPores();
	writeFrozenState();
	printSystemStats(frameNum);

	/*-- Calculate force once before starting the simulation. 
	 * This is done in order to avoid having to calculate the force 
	 * twice for each timestep */
	computeForce();
	cout << "\nStarting time integration. " << endl;
	cout << "-------------------------------------\n\n";
	clock_t start = clock(); clock_t end;
	int counter = 1;
	for(double t = dt; t <= finalT; t = counter*dt){
		doVerletIntegration();
		computeTemperature(false); 
		applyBerendsenThermostat();

		if(counter % dumpRate == 0){
			// State calculations
			frameNum++;
			computeEnergy();
			computeTemperature(true); 
			computePressure();
			end = clock();
			// Write out
			writeState(frameNum);
			cout << "Dumping .xyz at step " 
				<< counter << " of " << finalT/dt << endl;
			cout << "T: " << T << ", dt = " << dt 
				<< " , steps: " << finalT/dt << endl;
			cout << "Time since last dump: " 
				<< double(end - start)/CLOCKS_PER_SEC << endl;
			cout << "-------------------------------------" 
				<< endl << endl;
			// Start the clock again
			start = end;
		}
		counter++;
	}
	return 0;
}

void createCylindricPores(){
	double 	centerXY[2] = {0,0};
	double 	r = 0;

	for (int q = 0; q < natoms; q++) {
		vpFrozenAtoms[q] = true;
	}

	using namespace PeriodicBounds;
	for (int pore = 0; pore < poreCylCnt; pore++) {
		centerXY[0] = L*Random::ran0(seed);
		centerXY[1] = L*Random::ran0(seed);
		for (int i = 0; i < natoms; i++) {
			// Functions from PeriodicBounds
			rij[0] = getClosestDist(mpState[i][0],centerXY[0], L);
			rij[1] = getClosestDist(mpState[i][1],centerXY[1], L);

			r = sqrt(rij[0]*rij[0]+rij[1]*rij[1]);
			if (r < poreCylRadius) {
				vpFrozenAtoms[i] = false;
			}
		}
	}

	for (int j = 0; j < natoms; j++) {
		if (vpFrozenAtoms[j] == true) {
			mpState[j][3] = 0;
			mpState[j][4] = 0;
			mpState[j][5] = 0;
		}
	}
	cout << "Pore radius = " << poreRadius << " sigma." << endl;

}

void initLocalParam(){
	// Init simulation variables
	L 	= b*Nc;
	boxSize = 3; // 3 sigma in real units
	boxes 	= ceil(L/boxSize);
	mpA 	= matrix(natoms,3);
	//-- Divide the domain into boxes of size rc
	cBox = Cube(boxes,boxes,boxes); 
	vpLinkedList = new int[natoms];
	refreshBoxes();
	// Variables for pores
	vpFrozenAtoms	= new bool[natoms];
	if (usePores == false && useCylinders == false) {
		for (int q = 0; q < natoms; q++) {
			// Unfreeze all atoms
			vpFrozenAtoms[q] = false;
		}
	}

}


void printSystemStats(int frameNum){
	sprintf(tmpstr,"Starting point: %04d.xyz",frameNum);
	cout << tmpstr << endl; 

	cout << "\nL: " << L << ", boxSize: " << boxSize 
		<< ", boxes: " << boxes << endl;
	cout << "Density: " 
	     << mass*natoms / 
	     ((b*Nc*sigma*1E-10)*(b*Nc*sigma*1E-10)*(b*Nc*sigma*1E-10)) 
	     << " kg/m^3" << endl;
	cout << "Actual system length: " << (b*sigma*1E-10) << " m" << endl;
	cout << "Actual system volume: " 
	     << (b*sigma*1E-10)*(b*sigma*1E-10)*(b*sigma*1E-10) << " m^3"
	     << endl << endl;
}

void writeFrozenState(){
	/* This function writes the current state of the 
	 * system to a .xyz file in ASCII format. */
	char fileName[256];
    
	sprintf(fileName, "%sStates/Frozen.xyz", outputPath);
	//cout << fileName << endl;
	FILE *outFile;
    	int fatoms = 0;
	for(int i = 0; i < natoms; i++){
		if(vpFrozenAtoms[i]) fatoms++;
	}
    
	outFile = fopen(fileName, "w");
	fprintf(outFile, "%i\n", fatoms);
	fprintf(outFile, "Argon atom system state\n");
	for(int i = 0; i < natoms; i++){
		if(vpFrozenAtoms[i]){
		    //cout << "Found frozen atom" << endl;
		    fprintf(outFile, 
					"Ar %e %e %e %e %e %e\n", 
					mpState[i][0], 
					mpState[i][1], 
					mpState[i][2], 
					mpState[i][3], 
					mpState[i][4], 
					mpState[i][5]);
		}
	}
    
	fclose(outFile);
}

void createPores(){
	double 	center[3] = {0,0,0};
	double 	r 	= 0;

	for (int q = 0; q < natoms; q++) {
		vpFrozenAtoms[q] = false;
	}

	using namespace PeriodicBounds;
	for (int pore = 0; pore < poreCnt; pore++) {
		center[0] = L*Random::ran0(seed);
		center[1] = L*Random::ran0(seed);
		center[2] = L*Random::ran0(seed);
		for (int i = 0; i < natoms; i++) {
			// Functions from PeriodicBounds
			rij[0] = getClosestDist(mpState[i][0],center[0], L);
			rij[1] = getClosestDist(mpState[i][1],center[1], L);
			rij[2] = getClosestDist(mpState[i][2],center[2], L);

			r = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
			if (r < poreRadius) {
				vpFrozenAtoms[i] = true;
			}
		}
	}

	for (int j = 0; j < natoms; j++) {
		if (vpFrozenAtoms[j] == true) {
			mpState[j][3] = 0;
			mpState[j][4] = 0;
			mpState[j][5] = 0;
		}
	}
	cout << "Pore radius = " << poreRadius << " sigma." << endl;
}

void refreshBoxes(){
	//-- Refresh boxes and linked list;
	cBox = -1;
	int ix, iy, iz;
	for(int i = 0; i < natoms; i++){
		vpLinkedList[i] = -1;
		ix = int(mpState[i][0]/boxSize);
		iy = int(mpState[i][1]/boxSize);
		iz = int(mpState[i][2]/boxSize);
		vpLinkedList[i] = cBox(ix, iy, iz);
		cBox(ix, iy, iz) = i;  
	}
}

void writeState(int frameNum){
	/* This function writes the current state of the 
	 * system to a .xyz file in ASCII format. */
	char fileName[256];
	sprintf(fileName, "%sStates/%04d.xyz", outputPath ,frameNum);
	//cout << fileName << endl;
	FILE *outFile;

	outFile = fopen(fileName, "w");
	fprintf(outFile, "%i\n", natoms);
	fprintf(outFile, "Argon atom system state\n");
	for(int i = 0; i < natoms; i++){
		if(vpFrozenAtoms[i]){
		    //cout << "Found frozen atom" << endl;
		    fprintf(outFile, 
					"Af %e %e %e %e %e %e\n", 
					mpState[i][0], 
					mpState[i][1], 
					mpState[i][2], 
					mpState[i][3], 
					mpState[i][4], 
					mpState[i][5]);
		}
		else{
		    fprintf(outFile, 
					"Ar %e %e %e %e %e %e\n", 
					mpState[i][0], 
					mpState[i][1], 
					mpState[i][2], 
					mpState[i][3], 
					mpState[i][4], 
					mpState[i][5]);
		}
	}
    
	fclose(outFile);
}


void loadState(int frameNum){
	char fileName[256];
	sprintf(fileName, "%sStates/%04d.xyz", outputPath, frameNum);
	FILE* inFile;
	long lSize;
	inFile = fopen(fileName, "r");
	fseek(inFile, 0, SEEK_END);
	lSize = ftell(inFile);
	rewind(inFile);

	// Set the value of natoms
	fscanf(inFile, "%i", &natoms);

	// Create the state array
	mState = Matrix(natoms,6);
	mpState = mState.getArrayPointer();

	cout << "-----\n" << "Loading state " 
	    << fileName << endl 
	    << "Natoms: " << natoms << endl;
	char line[512];
	fgets(line, 512, inFile);
	fgets(line, 512, inFile);
	cout << line;
	// Read the position and velocity of the particles form the file
	for(int i = 0; i < natoms; i++){
		fscanf(inFile, "%s", line);
		fscanf(inFile, "%lf", &mpState[i][0]);
		fscanf(inFile, "%lf", &mpState[i][1]);
		fscanf(inFile, "%lf", &mpState[i][2]);
		fscanf(inFile, "%lf", &mpState[i][3]);
		fscanf(inFile, "%lf", &mpState[i][4]);
		fscanf(inFile, "%lf", &mpState[i][5]);
	}
	cout <<  "Done loading state\n-----\n";
    	fclose(inFile);
}

void applyForce(int i, int j){
	rij[0] = mpState[i][0] - mpState[j][0];
	rij[1] = mpState[i][1] - mpState[j][1];
	rij[2] = mpState[i][2] - mpState[j][2];
	double lHalf = L/2.0;

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
    
	/* Efficient calculation of the three needed numbers */
	r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	r2i = 1.0/r2;
	r6i = r2i * r2i * r2i;
	r12i = r6i * r6i;
    
	for(int k = 0; k < 3; k++){
		aij[k] = 24 * (2*r12i - r6i) * r2i * rij[k];
	}
	mpA[i][0] += aij[0];
	mpA[i][1] += aij[1];
	mpA[i][2] += aij[2];

	sumPressure += rij[0]*mpA[i][0] + rij[1]*mpA[i][1] + rij[2]*mpA[i][2];
}

void computeForce(){
	// The pressure calculation needs to be reset every force calculation
	sumPressure = 0; 
	refreshBoxes();
	int atom1 = 0; int atom2 = 0; 
	for (int i = 0; i < natoms; i++){
		mpA[i][0] = 0;
		mpA[i][1] = 0;
		mpA[i][2] = 0;
	}
    
    // Loop through all boxes
    for(int i = 0; i < boxes; i++){
        for(int j = 0; j < boxes; j++){
            for(int k = 0; k < boxes; k++){
                
                atom1 = cBox(i,j,k);

		// Skip frozen atoms
		if(vpFrozenAtoms[atom1] == true){
			atom1 = vpLinkedList[atom1];
		}
                
                while(atom1 != -1){
                    // Loop through all neighbouring boxes
                    for(int x = i-1; x <= i+1; x++){
                        for(int y = j-1; y <= j+1; y++){
                            for(int z = k-1; z <= k+1; z++){
                                
                                atom2 = cBox(x,y,z); // Allows periodic indicies
                                while(atom2 != -1){
                                    if(atom2 != atom1){
                                        applyForce(atom1,atom2);
                                    }
                               	    atom2 = vpLinkedList[atom2];
                                }
                                
                            }
                        }
                    }
	            atom1 = vpLinkedList[atom1];
                }
            }
        }
    }
}



void doVerletIntegration(){
	/* This function performs the time integration using the Velocity Verlet
	integration integration scheme. It takes a state array as argument which
	should contain a number of arrays on the form [x, y, z, vx, vy, vz],the
	number of atoms and the time-step*/
    
	for(int i = 0; i < natoms; i++){
		// Skip atom if frozen
		if (vpFrozenAtoms[i] == true) continue;

		mpState[i][3] += dt*0.5*mpA[i][0]; //vhalf
		mpState[i][4] += dt*0.5*mpA[i][1];
		mpState[i][5] += dt*0.5*mpA[i][2];

		mpState[i][0] += mpState[i][3]*dt;
		mpState[i][1] += mpState[i][4]*dt;
		mpState[i][2] += mpState[i][5]*dt;

		/* Insert periodic boundaries */
		PeriodicBounds::correctPos(mpState[i][0],L);
		PeriodicBounds::correctPos(mpState[i][1],L);
		PeriodicBounds::correctPos(mpState[i][2],L);
	}

	/* At this point, recalculate forces */
	computeForce();

	for(int i = 0; i < natoms; i++){ 
		// Skip atom if frozen
		if (vpFrozenAtoms[i] == true) continue;
		// Then set new velocities
		mpState[i][3] += dt*0.5*mpA[i][0];
		mpState[i][4] += dt*0.5*mpA[i][1];
		mpState[i][5] += dt*0.5*mpA[i][2];
        
	}
	/* Done with one step of the integration */
}

void computeEnergy(){
	/* This function computes the total sum of kinetic and potential
	energy. The kinetic energy in non-dimensional units is given by
	Ek = 0.5 * v * v. The kinetic energy is given by
	Ep = 4([1/r]^12 - [1/r]^6) */
	double E = 0; double Ek = 0; double Ep = 0;
	for(int i = 0; i < natoms; i++){
		Ep = 0;
		// Add the kinetic energy for the atom
		Ek = 0.5*(mpState[i][3]*mpState[i][3] 
				+ mpState[i][4]*mpState[i][4] 
				+ mpState[i][5]*mpState[i][5]);
		// Then compute the total potential energy for this atom;
		for(int j = 0; j < natoms && i!=j; j++){
			rij[0] = mpState[i][0] - mpState[j][0];
			rij[1] = mpState[i][1] - mpState[j][1];
			rij[2] = mpState[i][2] - mpState[j][2];

			/* Efficient calculation of the three needed numbers */
			r2 = rij[0]*rij[0] 
				+ rij[1]*rij[1] 
				+ rij[2]*rij[2];
			r2i = 1.0/r2;
			r6i = r2i * r2i * r2i;
			r12i = r6i * r6i;

			Ep += 4 * (r12i - r6i);
		}
		E += Ek + Ep;
	}
	cout << "Total energy: " << E << endl;

	char fileName[64];
	sprintf(fileName, "%s/Measurements/Energy.dat", outputPath);
	FILE *outFile;
	outFile = fopen(fileName, "a");
	fprintf(outFile, "%e\n", E);
	fclose(outFile);
}

void computeTemperature(bool write){
	// It is important the temperature is calculated before
	// pressure and the thermostat
	double sumTemp = 0;
	for(int i = 0; i < natoms; i++){
		sumTemp += 0.5*(mpState[i][3]*mpState[i][3] 
				+ mpState[i][4]*mpState[i][4] 
				+ mpState[i][5]*mpState[i][5]);
	}

	double newT = 2 * sumTemp / (3 * natoms * K_B);
	if (write){
		cout << "Temperature: " << newT*e0 << " Kelvin." << endl;
		char fileName[64];
		sprintf(fileName, "%s/Measurements/Temperature.dat", outputPath);
		FILE *outFile;
		outFile = fopen(fileName, "a");
		fprintf(outFile, "%e\n", newT*e0);
		fclose(outFile);
	}
	curT = newT*e0;
}

void computeRadialDistr(double lHalf){
    
    /* Use distance bins with width L/100 */
    int binCount[100];
    // Set all bincount to zero initially
    for(int i = 0; i < 100; i++){ binCount[i] = 0; }
    
    
    double binWidth = lHalf/100.;
    /* First calculate the distance between the atoms */
    double r;
    int bindex;
    for(int i = 0; i < natoms; i++){
        rij[0] = mpState[0][0] - mpState[i][0];
        rij[1] = mpState[0][1] - mpState[i][1];
        rij[2] = mpState[0][2] - mpState[i][2];
        double lHalf = L/2.0;
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
        
        r = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        if(r < lHalf){
            bindex = int(r/binWidth);
            //cout << "bindex: " << bindex << endl;
            binCount[bindex] += 1;
        }
    }
    
    char fileName[128];
    sprintf(fileName, "%s/Measurements/RadialDistr.dat", outputPath);
    FILE *outFile;
    outFile = fopen(fileName, "a");
    for(int i = 0; i < 100; i++){
        fprintf(outFile, "%i ", binCount[i]);
    }
    fprintf(outFile, "\n");
    fclose(outFile);
}

void computePressure(){
	// Calculate pressure
	sumPressure = sumPressure/2;
	double V = L*L*L;
	sumPressure = sumPressure/(3*V);
	sumPressure = sumPressure*mass*mass/(3*b*Nc*Nc*Nc*sigma*e0);

	double rhoKT = 4/(b*b*b)*K_B*curT;
	rhoKT = rhoKT/(sigma*sigma*sigma)*1e30;

	double pressure =  rhoKT + sumPressure;

	cout << "Pressure : " << pressure/1000 << " kPa." << endl;

	char fileName[128];
	sprintf(fileName, "%s/Measurements/Pressure.dat", outputPath);
	FILE *outFile;
	outFile = fopen(fileName, "a");
	fprintf(outFile, "%e ", pressure);
	fprintf(outFile, "\n");
	fclose(outFile);
}

void applyBerendsenThermostat(){
    double gamma = sqrt(1 + dt/tau*(targetT/curT - 1));
    //cout << "Gamma: " << gamma << endl;
    for(int i = 0; i < natoms; i++){
	// Skip if frozen
	if (vpFrozenAtoms[i] == true) continue;
        mpState[i][3] *= gamma;
        mpState[i][4] *= gamma;
        mpState[i][5] *= gamma;
    }
}

void applyAndersenThermostat(){
    double ran = 0;
    for(int i = 0; i < natoms; i++){
        if(vpFrozenAtoms[i] = true) continue;
        // Generate random number
        // if random number less than dt/tau -> reassign velocity
        ran = Random::ran2(seed);
        if(ran < dt/tau){
            ran = Random::gauss(seed);
            mpState[i][3] = ran*stdDev;
            ran = Random::gauss(seed);
            mpState[i][4] = ran*stdDev;
            ran = Random::gauss(seed);
            mpState[i][5] = ran*stdDev;
        }
    }   
}

