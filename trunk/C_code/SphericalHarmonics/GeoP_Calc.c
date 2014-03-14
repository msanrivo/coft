#include <math.h>
#include <stdlib.h>

// Physical Constants
#define  GM 3986004.415e+8
#define  R  6378136.3

void Legendre_V3(int N,long double NN, long double P[], long double Pn_1[], long double Pn_2[], long double x);
void initializeALF(long double P[], long double Pn_1[],long double Pn_2[],int N);

long double GeoP_Calc(long double r, long double phi, long double lambda, long double Cnm[], long double Snm[],int degree){

	// Declaring Required Variables for Algorithm

	/*
	 * n    = the degree of the ALF polynomial
	 *
	 * j    = the array location to read ALF value. It corresponds to the ORDER of ALF
	 *
	 * pos  = the array location for Mass Coefficients Cnm and Snm --> =  a counter
	 *
	 *
	 *U     = GeoPotential Value
	 *
	 *
	 * Recursions Sum is the value of the outermost sigma summation in the algorithm
	 *
	 * sumj is the value of the innermost sigma summation in the algorithm
	 *
	*/

	int n, j, pos;

	double U = 0;
	long double Recursions_Sum = 0;
	long double sumj = 0;

	//Initializing Mass Coefficients Position
	pos = 0;

	// Allocating Space for ALF Values

	/* P is the actual Legendre polynomial
	 * Pn_1 is the previous Legendre polynomial as required in the recursive algorithm
	 * Pn_2 is the 2 previous Legendre polynomial as required in the recursive algorithm
	 *
	 */

	long double P[degree];
	long double Pn_1[degree];
	long double Pn_2[degree];

	initializeALF(P,Pn_1,Pn_2,degree);


	for(n = 1; n <= degree; n++){

			Legendre_V3(n,(long double) n, P, Pn_1, Pn_2, sin(phi)); // Computing Legendre Polynomial of degree n and orders from 0 to n

			sumj = 0;

			for(j=0 ; j <= n; j++){

				pos = pos+1;
				sumj =   pow((R/r),n) * P[j] *( Cnm[pos] * cos(j*lambda) + Snm[pos] * sin(j*lambda) ) + sumj;

			}


		Recursions_Sum = Recursions_Sum + sumj;


		}

		U = (1 + Recursions_Sum) * GM/r;

		return U;

}


