#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Physical Constants

#define  GM 3986004.415e+8
#define  R  6378136.3
#define  pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#define  itggoce_array_size 29160

// File Functions
int Reading_Mass_Coefficients(long double C[], long double S[]);
long double GeoP_Calc(long double r, long double phi, long double lambda, long double Cnm[], long double Snm[],int degree);

int main(){

	printf("Initializing Program");

	// time variable needed to measure computational time

	clock_t begin, end;
	double time_spent;

	begin = clock(); // Start Measuring Time

	// Location Coordinates

	long double lambda = 0;
	long double phi = 0;
	long double r = 1.5*R;

	// GeoPotential Value

	long double U = 0;

	// Reading database for Mass Coefficients

	long double Cnm[itggoce_array_size];
	long double Snm[itggoce_array_size];
	int degree = 240; // Maximum Degree Acquired;


	Reading_Mass_Coefficients(Cnm,Snm);

	//Calling GeoPotential Calculator Function

	U = GeoP_Calc(r,phi,lambda,Cnm,Snm,degree);

	printf("\n------------------------------------\n");
	printf("\nGeoPotential at coordinates (r,lambda,phi) = (%.2f * R ,%.3f ,%.3f ) is U = %f\n",(double) r/R,(double) lambda, (double) phi, (double) U);
	printf("\n-------------------------------------\n");

	printf("\n\nFinished\n");

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	printf("Running Time %f",time_spent);



	return 0;
}
