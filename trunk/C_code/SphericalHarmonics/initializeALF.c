void initializeALF(long double P[], long double Pn_1[],long double Pn_2[],int N){


	int i = 0;

	for(i=0;i<=N;i++){

		P[i]    = 0;
		Pn_1[i] = 0;
		Pn_2[i] = 0;

	}

	// Algorithm Seed
	Pn_1[0] = 1;

}
