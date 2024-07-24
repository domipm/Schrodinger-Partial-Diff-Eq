#include<cmath>
#include<iostream>
#include"complex.h"

#define MAX 1000

#define PI 3.14159265358979323846
#define H 6.62607004e-34

int N, n_ciclos, n_temp;
double lam, k, s, norm = 0;
double V[MAX]; // Array of potential values

fcomplex fonda[MAX];                                        // Wave function (complex)
fcomplex A[MAX], alpha[MAX], B[MAX], beta[MAX], chi[MAX];   // Algorithm coefficients
fcomplex i = Complex(0,1);

using namespace std;

int main() {


    //  Output files
    FILE *out;
    out = fopen("fonda.dat", "w");  // Wavefunction modulus
    FILE *pot;
    pot = fopen("pot.dat", "w");    // Potential
    FILE *norma;
    norma = fopen("norm.dat", "w"); // Wavefunction norm


    //  Initial parameters

    N = 1000;       // X-axis divisions
    n_ciclos = 100; // Cycles to-do
    n_temp = 10000; // Number of temporal iterations
    lam = 0.8;      // Lambda parameter (height of potential)

    k = 2*PI*double(n_ciclos)/double(N);    // Rescaling factor K
    s = 1/(4*k*k);                          // Rescaling of temporal divisions

    //  Potential
    for (int j = 0; j < N; j++) {
        if ( (j < 3*N/5 ) && (j > 2*N/5) )  V[j] = lam*k*k; // Potential barrier of height LAMBDA*K*K
        else V[j] = 0.0;                                    // Null potential outside defined range
        fprintf(pot, "%.15lf\n", V[j]);                     // Write potential to file
    }
    //  Close file
    fclose(pot);

    //  Initial wave function
    for (int j = 0; j < N; j++) {
        if (j == 0 || j == (N)) fonda[j] = Complex(0,0);    // Boundary conditions (wavefunction zero at extremes)
        else {
            //  Obtain wave function at each position
            fonda[j] = Complex( cos(k*j) * exp( ( -8*double(4*j-N)*double(4*j-N) ) / double(N*N) ) , sin(k*j) * exp( ( -8*double(4*j-N)*double(4*j-N) ) / double(N*N) ) );
            norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));  // Norm of wave function
        }
    }

    //  Normalization
    for (int j = 0; j < N; j++) {
        fonda[j] = Cdiv(fonda[j], Complex(sqrt((norm)),0)); // Divide by norm
        fprintf(out, "%.15lf\n", Cabs(fonda[j]));           // Rewrite normalized wavefunction at each point
    }
    fprintf(out, "\n"); // Newline for new iterations
    

    // A_j^0, alpha vectors
    for (int j = 0; j < N; j++) A[j] = Complex( -2-V[j], 2/s );                                     // Vector A_j^0
    alpha[N-1] = Complex(0.0, 0.0);                                                                 // Boundary condition alpha_N-1 = 0
    for (int j = N-2; j > 0; j--) alpha[j-1] = Cdiv( Complex(-1.0, 0.0), Cadd(A[j], alpha[j] ) );   // Vector alpha_j


    //  Algorithm to obtain time evolution of wavefunction
    for (int t = 0; t < n_temp; t++) {

        norm = 0; // Zero norm at beginning of each iteration

        //  B_j, beta vectors
        for (int j = 0; j < N; j++) B[j] = Cdiv( Cmul( Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0) ); // Vector B_j
        beta[N-1] = Complex(0.0, 0.0);                                                                  // Boundary condition beta_N-1 = 0
        for (int j = N-2; j > 0; j--) beta[j-1] = Cdiv( Csub(B[j], beta[j]), Cadd(A[j], alpha[j]) );    // Vector beta_j

        //  Chi vector
        for (int j = 0; j < N-1; j++) chi[j+1] = Cadd( Cmul(alpha[j], chi[j]), beta[j] );

        //  Wave function
        for (int j = 0; j < N; j++) {
            fonda[j] = Csub( chi[j], fonda[j] );
            norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
        }

        fprintf(norma, "%d\t%.15lf\n", t, sqrt(norm));                          // Write square root of norm squared for each iteration
        for (int j = 0; j < N; j++) fprintf(out, "%.15lf\n", Cabs(fonda[j]));   // Write new modulus of wave function
        fprintf(out, "\n");                                                     // Newlines to separate iterations

    }

    //  Close files
    fclose(out);
    fclose(norma);

    return 0;

}