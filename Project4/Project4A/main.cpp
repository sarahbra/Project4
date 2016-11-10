#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include "lib.h"

using namespace std;
using namespace arma;

ofstream ofile;

inline int periodic(int i, int limit,int add) {
    return (i+limit+add) % (limit);
}

void initializeLattice(int Nspins, int** &SpinMatrix, int &Energy, int &MagneticMoment);
void Metropolis(int number_of_spins, long &idum, int &E, int &M, double *w, int **spin_matrix);
void output(int, int, double, double*);
double partition_function(double* w);
double numeric_heat_capacity(double* average, int MCcycles, double number_of_spins);
double analytical_heat_capacity(double number_of_spins);

int main()
{
    char *outfilename;
    long idum;
    int **spin_matrix, number_of_spins, mcs, E, M;
    double w[17], average[5], temperature, Z, Cv, Cv2;

    outfilename = "results.txt";
    ofile.open(outfilename);
    number_of_spins = 2;
    mcs = 100000;

    temperature = 1.0;
    spin_matrix = (int**) matrix(number_of_spins,number_of_spins,sizeof(int));

    idum = -1;
    E=M = 0;

    for (int de= -8; de <= 8; de++) w[de+8] = 0;
    for (int de= -8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    for(int i=0; i<5; i++) average[i] = 0;

    Z = partition_function(w);

    initializeLattice(number_of_spins, spin_matrix, E, M);

    for(int cycle=1; cycle<=mcs; cycle++) {
        Metropolis(number_of_spins,idum,E,M,w,spin_matrix);
        average[0] += E;
        average[1] += E*E;
        average[2] += M;
        average[3] += M*M;
        average[4] += fabs(M);
    }

    output(number_of_spins,mcs,temperature,average);
    ofile.close();
    return 0;
}

void initializeLattice(int Nspins, int** &SpinMatrix, int &Energy, int &MagneticMoment)
{
    for(int x =0; x<Nspins; x++){
        for(int y=0; y<Nspins; y++){
            SpinMatrix[x][y] = 1.0;
            MagneticMoment += (double) SpinMatrix[x][y];

        }
    }

    for(int x=0; x<Nspins; x++){
        for(int y=0; y<Nspins; y++) {
            Energy -= (double) SpinMatrix[x][y]*
                    (SpinMatrix[periodic(x,Nspins,-1)][y] +
                     SpinMatrix[x][periodic(y,Nspins, -1)]);
        }
    }
}


void Metropolis(int number_of_spins, long& idum, int& E, int& M, double *w, int **spin_matrix) {
    int i = 0;
    int t = 0;
    for(int x=0; x<number_of_spins; x++) {
        for(int y=0; y<number_of_spins; y++) {
            int ix = (int) (ran1(&idum)*(double)number_of_spins);
            int iy = (int) (ran1(&idum)*(double)number_of_spins);
            int dE = 2*spin_matrix[ix][iy]*(spin_matrix[ix][periodic(iy,number_of_spins,-1)] +
                    spin_matrix[periodic(ix,number_of_spins,-1)][iy] +
                    spin_matrix[ix][periodic(iy,number_of_spins,1)] +
                    spin_matrix[periodic(ix,number_of_spins,1)][iy]);
            if(ran1(&idum)<=w[dE+8]) {
                spin_matrix[ix][iy] *= -1;
                M += 2*spin_matrix[ix][iy];
                E += dE;
            }
        }
    }
}

double partition_function(double *w){
    double Z;

    Z = 0;
    for(int i = 0; i <17; i++){
        Z += w[i];
    }
    return Z;
}

double analytical_heat_capacity(double number_of_spins) {
    double Z, mean_E, C_v;
    // for J = beta = 1, the partition function reduces to
    Z = 4*cosh(8) + 12;
    mean_E = 32*sinh(8)/Z;
    C_v = ((256*cosh(8))/Z - mean_E*mean_E)*1/number_of_spins*1/number_of_spins;
    return C_v;
}

void output(int NSpins, int MCcycles, double temperature, double* ExpectationValues){
    double norm = 1.0/((double) (MCcycles));
    double E_ExpectationValues = ExpectationValues[0]*norm;
    double E2_ExpectationValues = ExpectationValues[1]*norm;
    double M_ExpectationValues = ExpectationValues[2]*norm;
    double M2_ExpectationValues = ExpectationValues[3]*norm;
    double Mabs_ExpectationValues = ExpectationValues[4]*norm;

    double C_v = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
    double susceptibility = (M2_ExpectationValues - M_ExpectationValues*M_ExpectationValues)/NSpins/NSpins;
    double M2_variance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;

    double C_v2 = analytical_heat_capacity(NSpins);

    cout << "Heat capacity numerical: " << C_v << ", heat capacity analytical: " << C_v2 << endl;
    cout << "Susceptibility numerical: " << susceptibility << ", susceptibility analytical: " << endl;

    ofile << setiosflags(ios::showpoint  |  ios::uppercase);
    ofile << setw(15) << setprecision(8) << NSpins;
    ofile << setw(15) << setprecision(8) << MCcycles;
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << E_ExpectationValues  ;
    ofile << setw(15) << setprecision(8) << E2_ExpectationValues ;
    ofile << setw(15) << setprecision(8) << M_ExpectationValues ;
    ofile << setw(15) << setprecision(8) << M2_ExpectationValues;
    ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues;
    ofile << setw(15) << setprecision(8) << ",";
}
