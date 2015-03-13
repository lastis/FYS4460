#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <random>
#include <cmath>
#include <time.h>
#include "Cube.cpp"
#include "Parameters.h"
using namespace std;

void Verlet_integration(double ***state, double ***a, int natoms, double dt, double L, Cube &box, int **ll);

void Write_state(double ***, int, int);

void Write_state_multiframe(double ***, int, int);

void Force(double ***a, double ***state, double L, int natoms, Cube &box, int **ll);

void Load_state(double ***state, int *natoms, int framenum);

void Init_lattice(double ***state, int Nc, double b);

void Init_velocities(double ***state, int natoms, double stddev);

void applyForce(double ***a, double ***state, double L);

void Force2(double ***a, double ***state, int natoms, Cube &box, int **ll);

void Compute_energy(double ***state, int natoms);

void Compute_temperature(double ***state, int natoms);

// State variables

extern const double     amu; 
extern const double     K_B; 
extern const double 	sigma;
extern const double 	b;
extern const double	    mass;
extern const double 	T;
extern const double	    e0;
extern const double	    v0;
extern const double	    stdDev;
// Simulation variables
extern const double 	dt;
extern const double 	finalT;
extern int 	            Nc;
extern int 	            atoms;
extern double 	        L;
extern double       	boxSize;
extern int 	            boxes;


void Refresh_boxes(Cube &box, int **ll, double ***state, int natoms, double rc){
    //-- Refresh boxes and linkedlist;
    box = -1;
    int ix, iy, iz;
    for(int i = 0; i < natoms; i++){
        (*ll)[i] = -1;
        ix = int((*state)[i][0]/rc);
        iy = int((*state)[i][1]/rc);
        iz = int((*state)[i][2]/rc);
        (*ll)[i] = box(ix, iy, iz);
        box(ix, iy, iz) = i;  
    }
}

int main(int nargs, char** argsv){
    
    //-- Get Lattice 
    int  framenum = atoi(argsv[1]);
    char tmpstr[64]; sprintf(tmpstr, "Starting point: %03d.xyz", framenum);
    double **pos;
    int natoms;
    cout << tmpstr << endl; 
    Load_state(&pos, &natoms, 0); 
    Nc = cbrt(natoms/4.);
    L = b*Nc;
    boxSize = L/2.;
    boxes = L/boxSize + 1;
    cout << "L: " << L << ", boxSize: " << boxSize << ", boxes: " << boxes << endl;
    //-- Constants
    /*double sigma = 3.405; //E-10
    double b = 5.260/sigma; // Aangstrom 
    double kb = 1.38E-23;
    double m = 39.948 * 1.66E-27;
    double T = 100.0;
    double eps0 = 119.8 * kb;
    double v0 = sqrt(eps0 / m); // 6.33E-3
    double stddev = sqrt(kb*T/m);
    double boxsize = Nc*b;
    
    cout << "b: " << b << endl;
    cout << "Velocity scale: " << v0 << endl;
    cout << "Standard deviation " << stddev << endl;
    */
   
    
    //-- Initialize array for calculating forces;
    double **a = new double*[natoms]; 
    for(int i = 0; i < natoms; i++){
        a[i] = new double[3]; 
        a[i][0] = 0; a[i][1] = 0; a[i][2] = 0;
    }
    
    //-- Divide the domain into boxes of size rc
    Cube box = Cube(boxes, boxes, boxes); box = -1;
    int *ll = new int[natoms];
    Refresh_boxes(box, &ll, &pos, natoms, boxSize);
    
    
    //-- Go through one list for checking
    bool tf = true;
    cout << box(0,1,0) << "  ";
    int p = box(0,0,0);
    while(tf){
        cout << ll[p] << "   ";
        p = ll[p];
        if(p==-1) tf=false;
    }
    
    
    /*-- Calculate force once before starting the simulation. This is done in 
    order to avoid having to calculate the force twice for each timestep */
    Force(&a, &pos, Nc*b, natoms, box, &ll);
    //Force2(&a, &pos, natoms, box, &ll);
    
    
    /* Perform time integration untill time T */
    cout << "Boxsize: " << boxSize << endl;
    int dumprate = 1;
    int counter = 0;
    cout << "T: " << T << ", dt = " << dt << " , steps: " << finalT/dt << endl;
    Compute_energy(&pos, natoms);
    Compute_temperature(&pos, natoms);
    for(double t = framenum*dumprate*dt; t<finalT; t+=dt){
        counter++;
        Verlet_integration(&pos, &a, natoms, dt, L, box, &ll);
        Compute_energy(&pos, natoms);
        Compute_temperature(&pos, natoms);
        if(counter % dumprate == 0){
            cout << "Dumping .xyz at step " << int(t/dt)+1 << " of " << finalT/dt << endl;
            framenum++;
            Write_state(&pos, natoms, framenum);
        }
    }
    
    
    
    /* To check stddev, mean and plot */
    //system("python plot_vel_distr.py");
    
    return 0;
}


void Write_state(double ***state, int natoms, int framenum){
    /* This function writes the current state of the system to a .xyz file in 
    ASCII format. */
    
    char filename[64];
    sprintf(filename, "./States/%03d.xyz", framenum);
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


void Write_state_multiframe (double ***state, int natoms, int framenum){
    /* This function writes the current state of the system to a .xyz file in 
    ASCII format. */
    
    char filename[64];
    sprintf(filename, "./States/MultiFrameState.xyz");
    //cout << filename << endl;
    if(framenum == 0){
        FILE *out = fopen(filename, "w");
        fclose(out);
    }
    FILE *outFile;
    outFile = fopen(filename, "a");
    fprintf(outFile, "%i\n", natoms);
    fprintf(outFile, "Argon atom system state\n");
    
    for(int i = 0; i < natoms; i++){
        fprintf(outFile, "Ar %e %e %e %e %e %e\n", (*state)[i][0], (*state)[i][1], (*state)[i][2], (*state)[i][3], (*state)[i][4], (*state)[i][5]);
    }
    fclose(outFile);
    
}

void Load_state(double ***state, int *natoms, int framenum){
    char filename[64];
    sprintf(filename, "./States/%03d.xyz", framenum);
    
    FILE *inFile;
    long lSize;
    inFile = fopen(filename, "r");
    fseek(inFile, 0, SEEK_END);
    lSize = ftell(inFile);
    rewind(inFile);
    
    fscanf(inFile, "%i", natoms);
    (*state) = new double*[(*natoms)];
    cout << "-----\n" << "Loading state " << filename << endl << "Natoms: " << (*natoms) << endl;
    char line[512];
    fgets(line, 512, inFile);
    fgets(line, 512, inFile);
    cout << line;
    for(int i = 0; i < (*natoms); i++){
        (*state)[i] = new double[6];
        fscanf(inFile, "%s", line);
        fscanf(inFile, "%lf", &(*state)[i][0]);
        fscanf(inFile, "%lf", &(*state)[i][1]);
        fscanf(inFile, "%lf", &(*state)[i][2]);
        fscanf(inFile, "%lf", &(*state)[i][3]);
        fscanf(inFile, "%lf", &(*state)[i][4]);
        fscanf(inFile, "%lf", &(*state)[i][5]);
    }
    cout << (*state)[0][0] << "  " << (*state)[0][1] << "  "  << (*state)[0][2] << "  "  << (*state)[0][3] << "  "  << (*state)[0][4] << "  "  << (*state)[0][5] << endl;
    cout << ".\n.\n.\n";
    
    cout <<  "Done loading state\n-----\n";
}

void applyForce(double ***a, double ***state, int i, int j, double L){
    double aij[3];
    double rij[3];
    rij[0] = (*state)[i][0] - (*state)[j][0];
    rij[1] = (*state)[i][1] - (*state)[j][1];
    rij[2] = (*state)[i][2] - (*state)[j][2];
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
    (*a)[i][0] += aij[0];
    (*a)[i][1] += aij[1];
    (*a)[i][2] += aij[2];
    (*a)[j][0] -= aij[0];
    (*a)[j][1] -= aij[1]; 
    (*a)[j][2] -= aij[2];    
}

void Force(double ***a, double ***state, double L, int natoms, Cube &box, int **ll){

    
    /* This function calculates the force between two given atoms */
    for (int i = 0; i < natoms; i++){
        (*a)[i][0] = 0;
        (*a)[i][1] = 0;
        (*a)[i][2] = 0;
    }
    for(int i = 0; i<natoms; i++){
        for(int j = i + 1; j < natoms; j++){
            /* Create the relative vector between the atoms */
            applyForce(a, state, i, j, L);
        }
    }
}

void Force2(double ***a, double ***state, int natoms, Cube &box, int **ll){
    //cout << "Starting force calculation" << endl;
    Refresh_boxes(box, ll, state, natoms, boxSize);
	int atom1 = 0; int atom2 = 0;
	int x = 0; int y = 0; int z = 0; 
    for (int i = 0; i < natoms; i++){
        (*a)[i][0] = 0;
        (*a)[i][1] = 0;
        (*a)[i][2] = 0;
    }
	for (int i = 0; i < boxes; i++) {
        //cout << "i: " << i << endl;
        for (int j = 0; j < boxes; j++) {
            //cout << "j: " << j << endl;
            for (int k = 0; k < boxes; k++) {
                //cout << "k: " << k << endl;
                // We choose box(i,j,k) and itterate through
                // the boxes around this box
                for (int x = i-1; x <= i+1; x++) {
                    //cout << "x: " << x << endl;
                    for (int y = j-1; y <= j+1; y++) {
                        //cout << "y: " << y << endl;
                        for (int z = k-1; z <= k+1; z++) {
                            //cout << "z: " << z << endl;
                            // This is the first index of the atom
                            // in the box(x,y,z)
                            atom1 = box(i,j,k);
                            atom2 = box(x,y,z);
                            //cout << "Atom1: " << atom1 << ", Atom2: " << atom2 << endl;
                            while (atom1 != -1) {
                                while (atom2 != -1) {
                                // Skip to next when the particle find
                                // itslf in the box
                                    //cout << "Atom1: " << atom1 << ", Atom2: " << atom2 << endl;
                                    if (atom1 == atom2){
                                        atom2 = (*ll)[atom2];
                                        continue;
                                    }	
                                    applyForce(a, state, atom1, atom2, L);
                                    // Done, now change j to the next 
                                    // particle in the box
                                    atom2 = (*ll)[atom2];
                                } // Second while loop ends here
                                // Done, now change j to the next 
                                // particle in the box
                                atom1 = (*ll)[atom1];
                                //cout << "Atom1: " << atom1 << endl;
                            }// First while loop ends here
                        } // z loop ends here
                    } // y loop ends here
                } // x loop ends here
            } // k loop ends here
        } // j loop ends here
	} // i loop ends here
    //cout << "Done with force calculation" << endl;
}



void Verlet_integration(double ***state, double ***a, int natoms, double dt, double L, Cube &box, int **ll){
    /* This function performs the time integration using the Velocity Verlet
    integration integration scheme. It takes a state array as argument which
    should contain a number of arrays on the form [x, y, z, vx, vy, vz],the
    number of atoms and the time-step*/
    
    
    
    for(int i = 0; i < natoms; i++){
        (*state)[i][3] += dt*0.5*(*a)[i][0]; //vhalf
        (*state)[i][4] += dt*0.5*(*a)[i][1];
        (*state)[i][5] += dt*0.5*(*a)[i][2];
        
        (*state)[i][0] += (*state)[i][3]*dt;
        (*state)[i][1] += (*state)[i][4]*dt;
        (*state)[i][2] += (*state)[i][5]*dt;
        
        /* Insert periodic boundaries */
        
        if( (*state)[i][0] > L){
            (*state)[i][0] -= L*(int((*state)[i][0]/L));
        }
        else if( (*state)[i][0] < 0){
            //cout << "Old pos: " << (*state)[i][0];
            (*state)[i][0] -= L*(int((*state)[i][0]/L)-1);
            //cout << " , new pos: " << (*state)[i][0] << endl;
        }
        if( (*state)[i][1] > L){
            (*state)[i][1] -= L*(int((*state)[i][1]/L));
        }
        else if( (*state)[i][1] < 0){
            (*state)[i][1] -= L*(int((*state)[i][1]/L)-1);
        }
        if( (*state)[i][2] > L){
            (*state)[i][2] -= L*(int((*state)[i][2]/L));
            //if((*state)[i][2] > L) cout << "Erronous position: z=" << (*state)[i][2] << endl;
        }
        else if( (*state)[i][2] < 0){
            (*state)[i][2] -= L*(int((*state)[i][2]/L)-1);
            //if((*state)[i][2] < 0) cout << "Erronous position: z=" << (*state)[i][2] << endl;
        }
    }
    /* At this point, recalculate forces */
    Force(a, state, L, natoms, box, ll);
    //Force2(a, state, natoms, box, ll);
    
    
    for(int i = 0; i < natoms; i++){ 
        // Then set new velocities
        (*state)[i][3] += dt*0.5*(*a)[i][0];
        (*state)[i][4] += dt*0.5*(*a)[i][1];
        (*state)[i][5] += dt*0.5*(*a)[i][2];
    }
    /* Done with one step of the integration */
}

void Compute_energy(double ***state, int natoms){
   /* This function computes the total sum of kinetic and potential
    energy. The kinetic energy in non-dimensional units is given by
    Ek = 0.5 * v * v. The kinetic energy is given by
    Ep = 4([1/r]^12 - [1/r]^6) */
    double rij[3];
    
    
    double E = 0; double Ek = 0; double Ep = 0;
    for(int i = 0; i < natoms; i++){
        Ep = 0;
        // Add the kinetic energy for the atom
        Ek = 0.5*((*state)[i][3]*(*state)[i][3] + (*state)[i][4]*(*state)[i][4] + (*state)[i][5]*(*state)[i][5]);
        // Then compute the total potential energy for this atom;
        for(int j = 0; j < natoms && i!=j; j++){
            rij[0] = (*state)[i][0] - (*state)[j][0];
            rij[1] = (*state)[i][1] - (*state)[j][1];
            rij[2] = (*state)[i][2] - (*state)[j][2];
            
            /* Efficient calculation of the three needed numbers */
            double r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
            double r2i = 1.0/r2;
            double r6i = r2i * r2i * r2i;
            double r12i = r6i * r6i;
            
            Ep += 4 * (r12i - r6i);
        }
        //cout << "Ek: " << Ek << ", Ep: " << Ep << endl;
        E += Ek + Ep;
    }
    cout << "Total energy: " << E << endl;
}

void Compute_temperature(double ***state, int natoms){
    
    double Ek = 0;
    for(int i = 0; i < natoms; i++){
        Ek += 0.5*((*state)[i][3]*(*state)[i][3] + (*state)[i][4]*(*state)[i][4] + (*state)[i][5]*(*state)[i][5]);
    }
    double newT = 2 * Ek / (3 * natoms * K_B);
    cout << "Temperature: " << newT*e0 << " Kelvin." << endl;
    
}
