#include <ctime>
#include <cmath>
#include <iostream>
using namespace std;

void fillRandomVec(double[], int);
void generateSolutions(double[], int);

int main()
{
    //Defining variables
    int n;
    cout << "Dimension of matrix:";
    cin >> n;
    double a[n-1], b[n], c[n-1], vecu[n], vecf[n];

    //Assigning random numbers to vector elements (other than zero, making the tridiagonal matrix unsolvable)
    srand((unsigned)time(NULL));
    fillRandomVec(b,n);
    fillRandomVec(a,n-1);
    fillRandomVec(c,n-1);
    generateSolutions(vecf,n);

    //Forwards substitution
    for (int k=1; k<=n; k++) {
        b[k] = b[k] - (a[k-1]*c[k-1])/b[k-1];
        vecf[k] = vecf[k] - (a[k-1]*vecf[k-1])/b[k-1];
    }

    //Backwards substitution
    vecu[n] = b[n];
    double temp;
    cout << "u" << endl;
    for (int k=n-1; k>=0; k--) {
        temp = c[k]*vecu[k+1];
        vecu[k] = 1.0/b[k]*((1.0/n*1.0/n*vecf[k]) - temp);
    }

    for (int k=0; k<n; k++) {
        cout << vecu[k] << endl;
    }
}

void fillRandomVec(double arr[], int n)  {
    for (int i=0; i<n; i++) {
        arr[i] = rand()%100 + 1;
    }
}

void generateSolutions(double arr[], int n) {
    double hstep = 1.0/n;
    double h = hstep;
    for (int i=0; i<n; i++) {
        arr[i] = 100*exp(-10*h);
        h += hstep;
        // cout << h << endl;
    }
}
