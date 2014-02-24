#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include "CPhys.h"
#include "Parameters.h"
using namespace std;
using namespace CPhys;

void verletIntegration(double**&   state, double**&   a, Cube& box, 
		        int*&     pointer, double&   sum);

void writeState(double**& state, int frameNum);

void loadState(Matrix& state, int frameNum);

void applyForce(double**& a, double**& state, int i, int j, double& sum);

void force(double**&       a, double**& state, Cube& box, 
	    int*&     pointer, double&     sum);

void computeEnergy(double**& state);

double computeTemperature(double**& state);

void computeRadialDistr(double**& state, double Lhalf);

void Compute_pressure(double sum, double curT);

void berendsenThermostat(double**& state, double T, double Tbath, double tau);

void refreshBoxes(Cube &box, int*& pointer, double**& state);

// State variables
extern const double 	sigma;
extern const double 	b;
extern const double	mass;
extern const double 	T;
extern const double	e0;
extern const double	v0;
extern const double	stdDev;
// Simulation variables
extern const double 	dt;
extern const double 	finalT;
extern const int	dumpRate;

int 	Nc;
int 	natoms;
double 	L;
double  boxSize;
int 	boxes;



int main(int nargs, char** argsv){
	/* Read starting point and initialize corresponding variables */
	int 	frameNum = atoi(argsv[1]);
	char 	tmpstr[64]; 
	Matrix	mState;
	double**	state;
	// Some state calculation variables
	double curT 	= 0;
	double sum 	= 0;
    
	loadState(mState,frameNum); 
	state 	= mState.getArrayPointer();

	// Print out
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
   
    
    
	//-- Initialize array for calculating forces;
	double **a = matrix(natoms,3);
    
	//-- Divide the domain into boxes of size rc
	Cube box = Cube(boxes,boxes,boxes); 
	box = -1; // Sets all entries to -1
	int* pointer = new int[natoms];
	refreshBoxes(box,pointer,state);
    
	/*-- Calculate force once before starting the simulation. 
	 * This is done in order to avoid having to calculate the force 
	 * twice for each timestep */
	force(a,state,box,pointer,sum);

    
	/* Perform time integration untill time T */
	int counter = 0;

	cout << "\nStarting time integration. " << endl;
	cout << "-------------------------------------\n\n";
	clock_t start = clock(); clock_t end;

	for(double t = dt; t<finalT; t+=dt){
		counter++;
		sum = 0;
		verletIntegration(state,a,box,pointer,sum);

		if(counter % dumpRate == 0){
			cout << "Dumping .xyz at step " 
				<< int(t/dt)+1 << " of " << finalT/dt << endl;
			frameNum++;
			writeState(state,frameNum);

			// State calculations

			cout << "Boxsize: " << boxSize << endl;
			cout << "T: " << T << ", dt = " << dt 
				<< " , steps: " << finalT/dt << endl;
			computeEnergy(state);
			curT = computeTemperature(state); 
			Compute_pressure(sum,curT);

			end = clock();
			cout << "Time since last dump: " 
				<< double(end - start)/CLOCKS_PER_SEC << endl;
			cout << "-------------------------------------" 
				<< endl << endl;
			start = end;
		}
	}
	return 0;
}

void refreshBoxes(Cube &box, int*& pointer, double**& state){
	//-- Refresh boxes and linked list;
	box = -1;
	int ix, iy, iz;
	for(int i = 0; i < natoms; i++){
		pointer[i] = -1;
		ix = int(state[i][0]/boxSize);
		iy = int(state[i][1]/boxSize);
		iz = int(state[i][2]/boxSize);
		pointer[i] = box(ix, iy, iz);
		box(ix, iy, iz) = i;  
	}
}

void writeState(double**& state, int frameNum){
	/* This function writes the current state of the 
	 * system to a .xyz file in ASCII format. */

	char fileName[64];
	sprintf(fileName, "../../res/States/%04d.xyz", frameNum);
	//cout << fileName << endl;
	FILE *outFile;

	outFile = fopen(fileName, "w");
	fprintf(outFile, "%i\n", natoms);
	fprintf(outFile, "Argon atom system state\n");
	for(int i = 0; i < natoms; i++){
		fprintf(outFile, 
				"Ar %e %e %e %e %e %e\n", 
				state[i][0], 
				state[i][1], 
				state[i][2], 
				state[i][3], 
				state[i][4], 
				state[i][5]);
	}
	fclose(outFile);
}


void loadState(Matrix& state, int frameNum){
	char fileName[64];
	sprintf(fileName, "../../res/States/%04d.xyz", frameNum);

	FILE* inFile;
	long lSize;
	inFile = fopen(fileName, "r");
	fseek(inFile, 0, SEEK_END);
	lSize = ftell(inFile);
	rewind(inFile);

	// Set the value of natoms
	fscanf(inFile, "%i", &natoms);

	state = Matrix(natoms,6);
	double** ppState = state.getArrayPointer();

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
		fscanf(inFile, "%lf", &ppState[i][0]);
		fscanf(inFile, "%lf", &ppState[i][1]);
		fscanf(inFile, "%lf", &ppState[i][2]);
		fscanf(inFile, "%lf", &ppState[i][3]);
		fscanf(inFile, "%lf", &ppState[i][4]);
		fscanf(inFile, "%lf", &ppState[i][5]);
	}

	cout << ppState[0][0] << "  " 
	     << ppState[0][1] << "  "  
	     << ppState[0][2] << "  "  
	     << ppState[0][3] << "  "  
	     << ppState[0][4] << "  "  
	     << ppState[0][5] << endl;
	cout << ".\n.\n.\n";
	cout <<  "Done loading state\n-----\n";

	// Init simulation variables
	Nc 	= cbrt(natoms/4.); 
	L 	= b*Nc;
	boxSize = 3; // 3 sigma in real units
	boxes 	= ceil(L/boxSize);
}

void applyForce(double**& a, double**& state, int i, int j, double& sum){
	double aij[3];
	double rij[3];
	rij[0] = state[i][0] - state[j][0];
	rij[1] = state[i][1] - state[j][1];
	rij[2] = state[i][2] - state[j][2];
	double Lhalf = L/2.0;

	if(rij[0] > Lhalf){
		rij[0] -= L;
	}
	else if(rij[0] < -Lhalf){
		rij[0] += L;
	}
	if(rij[1] > Lhalf){
		rij[1] -= L;
	}
	else if(rij[1] < -Lhalf){
		rij[1] += L;
	}
	if(rij[2] > Lhalf){
		rij[2] -= L;
	}
	else if(rij[2] < -Lhalf){
		rij[2] += L;
	}
    
	/* Efficient calculation of the three needed numbers */
	double r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	double r2i = 1.0/r2;
	double r6i = r2i * r2i * r2i;
	double r12i = r6i * r6i;
    
	for(int k = 0; k < 3; k++){
		aij[k] = 24 * (2*r12i - r6i) * r2i * rij[k];
	}
	a[i][0] += aij[0];
	a[i][1] += aij[1];
	a[i][2] += aij[2];
	//a[j][0] -= aij[0];
	//a[j][1] -= aij[1]; 
	//a[j][2] -= aij[2];

	sum += rij[0]*a[i][0] + rij[1]*a[i][1] + rij[2]*a[i][2];
}

void force(double**& a, double**& state, Cube &box, 
		int*& pointer, double &sum){

	refreshBoxes(box,pointer,state);
	int atom1 = 0; int atom2 = 0; 
	for (int i = 0; i < natoms; i++){
		a[i][0] = 0;
		a[i][1] = 0;
		a[i][2] = 0;
	}
    
    // Loop through all boxes
    for(int i = 0; i < boxes; i++){
        for(int j = 0; j < boxes; j++){
            for(int k = 0; k < boxes; k++){
                
                atom1 = box(i,j,k);
                
                while(atom1 != -1){
                    // Loop through all neighbouring boxes
                    for(int x = i-1; x <= i+1; x++){
                        for(int y = j-1; y <= j+1; y++){
                            for(int z = k-1; z <= k+1; z++){
                                
                                atom2 = box(x,y,z); // Allows periodic indicies
                                while(atom2 != -1){
                                    if(atom2 != atom1){
                                        applyForce(a,state,atom1,atom2,sum);
                                    }
                                    atom2 = pointer[atom2];
                                }
                                
                            }
                        }
                    }
                    atom1 = pointer[atom1];
                }
            }
        }
    }
}



void verletIntegration(double**& state, double**& a, Cube &box, 
		int*& pointer, double& sum){
	/* This function performs the time integration using the Velocity Verlet
	integration integration scheme. It takes a state array as argument which
	should contain a number of arrays on the form [x, y, z, vx, vy, vz],the
	number of atoms and the time-step*/
    
    
    
	for(int i = 0; i < natoms; i++){
		state[i][3] += dt*0.5*a[i][0]; //vhalf
		state[i][4] += dt*0.5*a[i][1];
		state[i][5] += dt*0.5*a[i][2];

		state[i][0] += state[i][3]*dt;
		state[i][1] += state[i][4]*dt;
		state[i][2] += state[i][5]*dt;

		/* Insert periodic boundaries */
		if( state[i][0] > L){
			state[i][0] -= L*(int(state[i][0]/L));
		}
		else if(state[i][0] < 0){
			state[i][0] -= L*(int(state[i][0]/L)-1);
		}
		if( state[i][1] > L){
			state[i][1] -= L*(int(state[i][1]/L));
		}
		else if(state[i][1] < 0){
			state[i][1] -= L*(int(state[i][1]/L)-1);
		}
		if( state[i][2] > L){
			state[i][2] -= L*(int(state[i][2]/L));
		}
		else if( state[i][2] < 0){
			state[i][2] -= L*(int(state[i][2]/L)-1);
		}
	}

	/* At this point, recalculate forces */
	force(a,state,box,pointer,sum);


	for(int i = 0; i < natoms; i++){ 
		// Then set new velocities
		state[i][3] += dt*0.5*a[i][0];
		state[i][4] += dt*0.5*a[i][1];
		state[i][5] += dt*0.5*a[i][2];
	}
	/* Done with one step of the integration */
}

void computeEnergy(double**& state){
	/* This function computes the total sum of kinetic and potential
	energy. The kinetic energy in non-dimensional units is given by
	Ek = 0.5 * v * v. The kinetic energy is given by
	Ep = 4([1/r]^12 - [1/r]^6) */
	double rij[3];
	double E = 0; double Ek = 0; double Ep = 0;
	for(int i = 0; i < natoms; i++){
		Ep = 0;
		// Add the kinetic energy for the atom
		Ek = 0.5*(state[i][3]*state[i][3] 
				+ state[i][4]*state[i][4] 
				+ state[i][5]*state[i][5]);
		// Then compute the total potential energy for this atom;
		for(int j = 0; j < natoms && i!=j; j++){
			rij[0] = state[i][0] - state[j][0];
			rij[1] = state[i][1] - state[j][1];
			rij[2] = state[i][2] - state[j][2];

			/* Efficient calculation of the three needed numbers */
			double r2 = rij[0]*rij[0] 
				+ rij[1]*rij[1] 
				+ rij[2]*rij[2];
			double r2i = 1.0/r2;
			double r6i = r2i * r2i * r2i;
			double r12i = r6i * r6i;

			Ep += 4 * (r12i - r6i);
		}
		E += Ek + Ep;
	}
	cout << "Total energy: " << E << endl;

	char fileName[64];
	sprintf(fileName, "../../res/Measurements/Energy.dat");
	FILE *outFile;
	outFile = fopen(fileName, "a");
	fprintf(outFile, "%e\n", E);
	fclose(outFile);
}

double computeTemperature(double**& state){
    
	double Ek = 0;
	for(int i = 0; i < natoms; i++){
		Ek += 0.5*(state[i][3]*state[i][3] 
				+ state[i][4]*state[i][4] 
				+ state[i][5]*state[i][5]);
	}
	double newT = 2 * Ek / (3 * natoms * K_B);
	cout << "Temperature: " << newT*e0 << " Kelvin." << endl;

	char fileName[64];
	sprintf(fileName, "../../res/Measurements/Temperature.dat");
	FILE *outFile;
	outFile = fopen(fileName, "a");
	fprintf(outFile, "%e\n", newT*e0);
	fclose(outFile);
	return newT*e0;
}

void computeRadialDistr(double**& state, double Lhalf){
    
    /* Use distance bins with width L/100 */
    int binCount[100];
    // Set all bincount to zero initially
    for(int i = 0; i < 100; i++){ binCount[i] = 0; }
    
    
    double binWidth = Lhalf/100.;
    /* First calculate the distance between the atoms */
    double rij[3];
    double r;
    int bindex;
    for(int i = 0; i < natoms; i++){
        rij[0] = state[0][0] - state[i][0];
        rij[1] = state[0][1] - state[i][1];
        rij[2] = state[0][2] - state[i][2];
        double Lhalf = L/2.0;
        if(rij[0] > Lhalf){
            rij[0] -= L;
        }
        else if(rij[0] < -Lhalf){
            rij[0] += L;
        }
        if(rij[1] > Lhalf){
            rij[1] -= L;
        }
        else if(rij[1] < -Lhalf){
            rij[1] += L;
        }
        if(rij[2] > Lhalf){
            rij[2] -= L;
        }
        else if(rij[2] < -Lhalf){
            rij[2] += L;
        }
        
        r = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        if(r < Lhalf){
            bindex = int(r/binWidth);
            //cout << "bindex: " << bindex << endl;
            binCount[bindex] += 1;
        }
    }
    
    char fileName[128];
    sprintf(fileName, "../../res/Measurements/RadialDistr.dat");
    FILE *outFile;
    outFile = fopen(fileName, "a");
    for(int i = 0; i < 100; i++){
        fprintf(outFile, "%i ", binCount[i]);
    }
    fprintf(outFile, "\n");
    fclose(outFile);
}

void Compute_pressure(double sum, double curT){
	// Calculate pressure
	sum = sum/2;
	double V = L*L*L;
	sum = sum/(3*V);
	sum = sum*mass*mass/(3*b*Nc*Nc*Nc*sigma*e0);

	double rhoKT = 4/(b*b*b)*K_B*curT;
	rhoKT = rhoKT/(sigma*sigma*sigma)*1e30;

	double pressure =  rhoKT + sum;

	cout << "Pressure : " << pressure/1000 << " kPa." << endl;

	char fileName[128];
	sprintf(fileName, "../../res/Measurements/Pressure.dat");
	FILE *outFile;
	outFile = fopen(fileName, "a");
	fprintf(outFile, "%e ", pressure);
	fprintf(outFile, "\n");
	fclose(outFile);
}

void berendsenThermostat(double**& state, double T, double Tbath, double tau){
    double gamma = sqrt(1 + dt/tau*(Tbath/T - 1));
    //cout << "Gamma: " << gamma << endl;
    for(int i = 0; i < natoms; i++){
        state[i][3] *= gamma;
        state[i][4] *= gamma;
        state[i][5] *= gamma;
    }
}

