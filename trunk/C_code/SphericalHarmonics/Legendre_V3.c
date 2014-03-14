#include <math.h>

void Legendre_V3(int N,long double NN, long double P[], long double Pn_1[], long double Pn_2[], long double x){

	// Defining loop variable i
	int i = 0;

	// Defining order of the polynomial
	int m[N];

	for(i = 0; i < N; i++){

		// Initializing

		m[i] = i;

	}

	// Defining chi_1 and chi_2 parameters

	long double chi_1 = 0;
	long double chi_2 = 0;


	// Computing Recursions

	if ( N == 1 ){

		chi_1 = sqrt( 3 );

		P[0] = 1/NN *  (2*NN-1) * x * chi_1 * Pn_1[0];
	}
	else {

		for(i = 0; i < N; i++){

			chi_1 = sqrt(         (2*NN + 1) * (NN - m[i])  /  (     (2*NN - 1) * (NN + m[i])    )         );
			chi_2 = sqrt(         (2*NN + 1) * (NN - m[i]) * (NN - m[i] -1) / (    (2*NN - 3) * ( NN + m[i]) * (NN + m[i] - 1)   )            );

			P[i] = ( 1/(NN-m[i]) ) * ( (2*NN-1) * x * Pn_1[i] * chi_1 - (NN + m[i] -1) * Pn_2[i] * chi_2 );

		}

	}


	P[N] = sqrt( (2*NN + 1) / (2*NN) ) * sqrt( 1 - x*x ) * Pn_1[N-1];

	// Copying array

	for(i = 0; i <= N ; i++){

		Pn_2[i] = Pn_1[i];
		Pn_1[i] = P[i];

	}




}

