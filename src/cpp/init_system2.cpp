#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include "CPhys.h"
#include "Parameters2.h"
using namespace std;
using namespace CPhys;


/////////////////////////////////////////////////
// 		Init State Variables
/////////////////////////////////////////////////

extern const double 	sigma;
extern const double 	b;
extern const double 	mass;
extern const double 	T;
extern const double 	e0;
extern const double 	v0;
extern const double 	stdDev;
extern const int	Nc;


void Init_velocities(double ***state, int natoms, double stddev){
    /* Create random-number generator */
	long seed = 123;
	double ran = 0;
    /* Now set random velocities for each atom from the normal distr */
    for(int i = 0; i < natoms; i++){
	ran = Random::gauss(seed);
        (*state)[i][3] = ran*stddev;
	ran = Random::gauss(seed);
        (*state)[i][4] = ran*stddev;
	ran = Random::gauss(seed);
        (*state)[i][5] = ran*stddev;
        //cout << distr(generator)*stddev << endl;
    }
    
    /* Now remove any linear momentum */
    double usum = 0; double vsum = 0; double wsum = 0;
    for(int i = 0; i < natoms; i++){
        usum += (*state)[i][3];
        vsum += (*state)[i][4];
        wsum += (*state)[i][5];
    }  
    double du = usum/natoms; double dv = vsum/natoms; double dw = wsum/natoms;
    for(int i = 0; i < natoms; i++){
        (*state)[i][3] -= du;
        (*state)[i][4] -= dv;
        (*state)[i][5] -= dw;
        
    }
    
}
void Init_lattice(double ***state, int Nc, double b){
    
    /* Loop through each lattice dimension */
    for(int i = 0; i < Nc; i++){
        for(int j = 0; j < Nc; j++){
            for(int k = 0; k < Nc; k++){
                int indx = i*Nc*Nc*4 + j*Nc*4 + k*4;
                /* Create vectors (6 indices) */
                (*state)[indx] = new double[6];
                (*state)[indx + 1] = new double[6];    
                (*state)[indx + 2] = new double[6];
                (*state)[indx + 3] = new double[6];
                    
                /* Then set the position of each atom */
                (*state)[indx + 0][0] = i*b;
                (*state)[indx + 0][1] = j*b;
                (*state)[indx + 0][2] = k*b;
            
                (*state)[indx + 1][0] = i*b + b/2.;
                (*state)[indx + 1][1] = j*b + b/2.;
                (*state)[indx + 1][2] = k*b;
            
                (*state)[indx + 2][0] = i*b;
                (*state)[indx + 2][1] = j*b + b/2.;
                (*state)[indx + 2][2] = k*b + b/2.;
                
                (*state)[indx + 3][0] = i*b + b/2.;
                (*state)[indx + 3][1] = j*b;
                (*state)[indx + 3][2] = k*b + b/2.;
            }
        }
    }
    cout << (*state)[Nc*Nc*Nc*4 - 1][0] << "   " << (*state)[Nc*Nc*Nc*4 - 1][1] << "   "  << (*state)[Nc*Nc*Nc*4 - 1][2] << endl;
}

void Write_state(double ***state, int natoms, int framenum){
    /* This function writes the current state of the system to a .xyz file in 
    ASCII format. */
    
    char filename[64];
    sprintf(filename, "../../res/States/%04d.xyz", framenum);
    //cout << filename << endl;
    FILE *outFile;
    
    outFile = fopen(filename, "w");
    fprintf(outFile, "%i\n", natoms);
    fprintf(outFile, "Argon atom system state\n");
    for(int i = 0; i < natoms; i++){
        fprintf(outFile, "Ar %e %e %e %e %e %e\n", (*state)[i][0], (*state)[i][1], (*state)[i][2], (*state)[i][3], (*state)[i][4], (*state)[i][5]);
    }
    fclose(outFile);
}


int main(int nargs, char** argsv){
    /* Check if user provided number of cells */
    int natoms = Nc*Nc*Nc*4;
    cout << "Specified " << Nc << " number of cells in each dimension" << endl;
    cout << "This gives "<< natoms << " atoms." << endl;
    cout << "Density: " << mass*natoms / 
	    ((b*Nc*sigma*1E-10)
	     *(b*Nc*sigma*1E-10)
	     *(b*Nc*sigma*1E-10)) 
	    << " kg/m^3" << endl;
    cout << "b: " << b << endl;
    double **pos = new double*[natoms];
    Init_lattice(&pos, Nc, b);
    Init_velocities(&pos, natoms, stdDev/v0);
    Write_state(&pos, natoms, 0);
    
}
