#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include "lib.h"
#include "mpi.h"


using namespace std;
using namespace arma;

ofstream ofile;
inline int periodic(int i, int limit,int add) {
    return (i+limit+add) % (limit);
}

void initializeLattice(int Nspins, long &idum, int** &SpinMatrix, double &Energy, double &MagneticMoment, int ordered);
void Metropolis(int number_of_spins, long &idum, double &E, double &M, double *w, int **spin_matrix, int& accepted_conf);
void output(int NSpins, int MCcycles, double temperature, double* ExpectationValues, int& accepted_conf);
double partition_function();

double numeric_heat_capacity(double* average, int MCcycles, int number_of_spins);
double heat_capacity(int number_of_spins);
double susceptibility(int number_of_spins);

int main(int argc, char* argv[])
{
    char outfilename[60];
    long idum;
    int **spin_matrix, number_of_spins, mcs, count, accepted_configurations;
    double w[17], average[5], temperature, E, M;
    int NProcesses, RankProcess;


    number_of_spins = 140;


    //int runs = 1;
    accepted_configurations = 0;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &NProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &RankProcess);


    mcs = 700000;

    sprintf(outfilename,"results4e_Spins%d_mcs%d.txt", number_of_spins,mcs);
    ofile.open(outfilename);

    ofile <<"NSpins" << "   Cycles   " <<"   Temperature  " << "E" << "Mabs "<< "  C_v" << "  chi" << endl;
    spin_matrix = (int**) matrix(number_of_spins,number_of_spins,sizeof(int));

    idum = -1;
    count = 0;


    int Energy[mcs-6000];
    int j =1;


    for(temperature=2.0; temperature <2.31; temperature += 0.01 ){
        for (int j=0;j<5;j++) {
           average[j] = 0.0;
        }
        cout << "temperature" << temperature <<"   " << RankProcess << endl;
        E=M=0;
        initializeLattice(number_of_spins, idum, spin_matrix, E, M, 0);
        for (int de= -8; de <= 8; de++) w[de+8] = 0;
        for (int de= -8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
        for(int i=0; i<5; i++) average[i] = 0;


        for(int cycle=1; cycle<=mcs; cycle++) {

                Metropolis(number_of_spins,idum,E,M,w,spin_matrix,accepted_configurations);
                    average[0] += E;
                    average[1] += E*E;
                    average[2] += M;
                    average[3] += M*M;
                    average[4] += fabs(M);


        }

        double TotalExpectation[5];
        for (int i=0; i<5; i++) TotalExpectation[i] = 0;

        for(int i=0; i<5; i++){
            MPI_Reduce(&average[i], &TotalExpectation [i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        if ( RankProcess == 0) output(number_of_spins,mcs,temperature,average, accepted_configurations);


    }

    ofile.close();

    MPI_Finalize ();
    return 0;

}

void initializeLattice(int Nspins, long &idum, int** &SpinMatrix, double &Energy, double &MagneticMoment, int ordered)
{
    for(int x =0; x<Nspins; x++){
        for(int y=0; y<Nspins; y++){
            if(ordered==0) {
                if(ran1(&idum)<=0.4999999999) {
                    SpinMatrix[x][y] = -1;
                } else {
                    SpinMatrix[x][y] = 1;
                }

            MagneticMoment += (double)SpinMatrix[x][y];

            } else {

            SpinMatrix[x][y] = 1.0;
            MagneticMoment += (double)SpinMatrix[x][y];
            }
        }
    }

    for(int x=0; x<Nspins; x++){
        for(int y=0; y<Nspins; y++) {
            Energy -= (double)SpinMatrix[x][y]*
                    (SpinMatrix[periodic(x,Nspins,-1)][y] +
                     SpinMatrix[x][periodic(y,Nspins, -1)]);
        }
    }
}


void Metropolis(int number_of_spins, long& idum, double& E, double& M, double *w, int **spin_matrix, int& accepted_conf) {

    for(int x=0; x<number_of_spins; x++) {
        for(int y=0; y<number_of_spins; y++) {
            int ix = (int) (ran1(&idum)*(double)number_of_spins);
            int iy = (int) (ran1(&idum)*(double)number_of_spins);
            int dE = 2*spin_matrix[ix][iy]*(spin_matrix[ix][periodic(iy,number_of_spins,-1)] +
                    spin_matrix[periodic(ix,number_of_spins,-1)][iy] +
                    spin_matrix[ix][periodic(iy,number_of_spins,1)] +
                    spin_matrix[periodic(ix,number_of_spins,1)][iy]);

            if(ran1(&idum)<=w[dE+8]) {

                spin_matrix[ix][iy] *= -1.0;
                M += (double)2*spin_matrix[ix][iy];
                E += (double)dE;

                accepted_conf += 1;

            }
       }
   }
}

double partition_function() {
    double Z = 4*cosh(8) + 12;
    return Z;
}

double heat_capacity(int number_of_spins) {
    double Z, E, C_v, E2;

    // for J = beta = 1, the partition function reduces to
    Z = partition_function();
    E = 32*sinh(8)/Z;
    E2 = (256*cosh(8))/Z;
    C_v = (E2 - E*E)*1/number_of_spins*1/number_of_spins;
    return C_v;
}

double susceptibility(int number_of_spins) {
    double chi, m, m2, abs_m, Z;
    Z = partition_function();
    m = 8*exp(8)+16;
    abs_m = fabs(m)/Z;
    m2 = (8*(exp(8)+1))/(cosh(8)+3);
    chi = (m2 - abs_m*abs_m)*1/number_of_spins*1/number_of_spins;
    return chi;
}

void output(int NSpins, int MCcycles, double temperature, double* ExpectationValues, int& accepted_conf){
    double norm = 1.0/((double) (MCcycles));
    double E_ExpectationValues = ExpectationValues[0]*norm;
    double E2_ExpectationValues = ExpectationValues[1]*norm;
    double M_ExpectationValues = ExpectationValues[2]*norm;
    double M2_ExpectationValues = ExpectationValues[3]*norm;
    double Mabs_ExpectationValues = ExpectationValues[4]*norm;
    double E_variance = (E2_ExpectationValues - E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
    double C_v = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins/temperature;

    double M_variance = (M2_ExpectationValues - M_ExpectationValues*M_ExpectationValues)/NSpins/NSpins;
    double chi = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins/temperature;

    double C_v2 = heat_capacity(NSpins);
    double chi2 = susceptibility(NSpins);
    E_ExpectationValues = E_ExpectationValues/NSpins/NSpins;
    Mabs_ExpectationValues = Mabs_ExpectationValues/NSpins/NSpins;

    //double C_v2 = heat_capacity(NSpins);
    //ouble chi2 = susceptibility(NSpins);

    //cout << "Heat capacity numerical: " << C_v << ", heat capacity analytical: " << C_v2 << endl;
    //cout << "Susceptibility numerical: " << chi << ", susceptibility analytical: " << chi2 << endl;


    ofile << setiosflags(ios::showpoint  |  ios::uppercase);
    ofile << NSpins<<",";
    ofile <<  MCcycles<< ",";
    ofile << temperature << ",";

    ofile << E_ExpectationValues <<",";
    //ofile << setw(15) << setprecision(8) << E2_ExpectationValues ;
    //ofile << setw(15) << setprecision(8) << M_ExpectationValues ;
    //ofile << setw(15) << setprecision(8) << M2_ExpectationValues;
    ofile <<Mabs_ExpectationValues<<",";
    //ofile << setw(20) << setprecision(8) << accepted_conf << endl;
    //ofile << setw(15) << setprecision(8) << E_variance << endl;
    ofile << C_v <<",";
    ofile << chi << endl;



}
