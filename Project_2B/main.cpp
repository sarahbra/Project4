#include<iostream>
#include <cmath>
using namespace std;

//Generating Schroedinger-matrix
void generate_matrix (double ** A, int n, double max_rho, double w_r) {
    double h = max_rho/(double) n;
    double hh = h*h;
    double rho = h;
    double nondiag = -1.0/(h*h);
    for (int i=0; i<(n-2); i++) {
        for (int j=0; j<(n-2); j++) {
            if (i==j) {
                if (w_r == 0) {
                    A[i][j] = 2.0/hh + rho*rho;
                } else {
                    A[i][j] = 2.0/hh + w_r*w_r*rho*rho + 1.0/rho;
                }
            } else if (j==(i+1) || j==(i-1)) {
                A[i][j] = nondiag;
            } else {
                A[i][j] = 0.0;
            }
        }
        rho += h;
        cout << "rho = " << rho << endl;
    }
    return;
}

void rotate (double ** A, double ** R, int k, int l, int n) {
    double c, s;
    if(A[k][l] != 0) {
        double t, tau;
        tau = (A[l][l]-A[k][k])/(2.0*A[k][l]);
        if (tau > 0) {
            t = 1.0/(tau + sqrt(1.0+tau*tau));
        } else {
            t = -1.0/(-tau + sqrt(1.0+tau*tau));
        }
        c = 1.0/sqrt(1+t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }
    double a_ll, a_kk, a_il, a_ik, r_il, r_ik;
    a_ll = A[l][l];
    a_kk = A[k][k];

    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][k] = 0.0;
    A[k][l] = 0.0;

    for ( int i = 0; i < n; i++ ) {
        //cout << "A" << i << i << " = " << A[i][i] << endl;
        if ( i != k && i != l ) {
            a_il = A[i][l];
            a_ik = A[i][k];
            A[i][l] = c*a_ik + s*a_il;
            A[l][i] = A[i][l];
            A[i][k] = c*a_il - s*a_ik;
            A[k][i] = A[i][k];
        }
        r_il = R[i][l];
        r_ik = R[i][k];
        R[i][l] = c*r_il + s*r_ik;
        //cout << "R" << i << l  << " = " << R[i][l] << endl;
        R[i][k] = c*r_ik - s*r_il;
    }
    return;
}

double maxnondiag (double ** A, int * k, int * l, int n) {
    double max = 0.0;
    double temp;
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++){
            temp = sqrt(A[i][j]*A[i][j]);
            // cout << "temp = " << temp << endl;
            if (temp > max) {
                max = temp;
                *l = i;
                // cout << "l = " << *l << endl;
                *k = j;
                // cout << "k = " << *k << endl;
            }
        }
    }
    return max;
}

void jacobi_method (double ** A, double ** R, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                R[i][j] = 1.0;
            } else {
                R[i][j] = 0.0;
            }
        }
    }

    int l, k;
    double epsilon = 1.0e-8;
    int iterations = 0;
    int max_iterations = n*n*n;
    double max_nondiag = maxnondiag(A, &k, &l, n);
    double eigenvalues[n];
    // cout << max_nondiag << endl;

    while ((iterations < max_iterations) && (sqrt(max_nondiag*max_nondiag) > epsilon)) {
        max_nondiag = maxnondiag(A, &k, &l, n);
        rotate (A, R, k, l, n);
        iterations++;
    }
    cout << "Number of iterations: " << iterations << endl;

    //finding the eigenvalues
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = A[i][i];
        cout << eigenvalues[i] << endl;
    }

    return;
}


//void fill_matrix(double ** A, int n)  {
//    for (int i=0; i<n; i++) {
//        for (int j=0; j<n; j++) {
//            if (i==(j-1)) {
//                A[i][j] = 10.0;
//            } else {
//                A[i][j] = rand()%9 + 1;
//            }
//        }
//    }
//    return;
//}


int main() {

    double ** A;
    int n = 350;
    double ** R;
    A = new double*[n-2];
    R = new double*[n-2];
        for (int i = 0; i<(n-2); ++i) {
            A[i] = new double[n-2];
            R[i] = new double[n-2];
        }
    int rho = 10;

    generate_matrix(A,n,rho,0);

    jacobi_method(A,R,(n-2));

    //int k,l;
    //double ** A2;
    //A2 = new double*[n];
    //for (int i = 0; i<n; i++) {
    //    A2[i] = new double[n];
    //}

    //fill_matrix(A,n);
    //double max_nondiag2 = maxnondiag(A,&k,&l,n);
    //cout << max_nondiag2 << endl;

}
