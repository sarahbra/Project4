#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <fstream>

using namespace std;
using namespace arma;

ofstream ofile;

inline int periodic(int i, int limit,int add) {
    return (i+limit+add) % (limit);
}

void initializeLattice(int Nspins, mat &SpinMatrix, double& Energy, double MagneticMoment);
void Metropolis(int number_of_spins, long &idum, int &E, int &M, double *w, int **spin_matrix);
void output(int, int, double, vec);

int main()
{
    return 0;
}

void initializeLattice(int Nspins, mat &SpinMatrix, int &Energy, int &MagneticMoment);
{
    for(int x =0; x<Nspins; x++){
        for(int y=0; y<Nspins; y++){
            SpinMatrix(x,y) = 1.0;
            MagneticMoment += (double) SpinMatrix(x,y);

        }
    }
    for(int x=0; x<Nspins; x++){
        for(int y=o; y<Nspins; y++)
            Energy -= (double) SpinMatrix(x,y)*
                    (SpinMatrix(periodic(x,Nspins,-1),y) +
                     SpinMatrix(x,periodic((y,Nspins, -1)));

}
}


void Metropolis(int number_of_spins, long& idum, int& E, int& M, double *w, int **spin_matrix) {
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
void output(int NSpins, int MCcycles, double temperature, double ExpectationValues){
    double norm = 1.0/((double) (MCcycles));
    double E_ExpectationValues = ExpectationValues(0)*norm;
    double E2_ExpectationValues = ExpectationValues(1)*norm;
    double M_expectationValues = ExpectationValues(2)*norm;
    double M2_ExpectationValues = ExpectationValues(3)*norm;
    double Mavs_Expectationvalues = ExpectationValues(4)*norm;

    //Variansen

    ofile << setiosflags(ios::showpoint  |  ios::uppercase);
    ofile << setw(15) << setprecision(8) <<  ;
    ofile << setw(15) << setprecision(8) << ;
    ofile << setw(15) << setprecision(8) << ;
    ofile << setw(15) << setprecision(8) << ;
    ofile << setw(15) << setprecision(8) << ;
    ofile << setw(15) << setprecision(8) << ;

}
