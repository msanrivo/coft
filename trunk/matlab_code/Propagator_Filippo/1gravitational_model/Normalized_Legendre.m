function [P, Pn_1, Pn_2]  = Normalized_Legendre(N,Pn_1,Pn_2,sinphi,cosphi)

%{
    Recursive Computation Of The Associated Legendre Polynomials of Order N

    Author: Diego García Pardo (UNIVERSITY CARLOS III OF MADRID)
    
    
    Requires as Input lower order Legendre Polynomials P @ N-1 and P @ N-2
    
    Also requires the introduction of the values of the sine and cosine of
    the longitude angle (measured from the x-y plane (horizontal with
    z-axis over the north pole)
    
    Therefore this function is ONLY VALID for greater or equal than 2 polynomial

%}

P = nan(N+1,1);
m = (0:N-2)'; 

chi_1 = sqrt(  (2*N+1).*(N-m) ./( (2*N-1).*(N+m) )  );
chi_2 = sqrt(  (2*N+1).*(N-m).*(N-m-1) ./ ( (2*N-3).*(N+m).*(N+m-1) )  );

chi_3 = sqrt(2*N+1);
chi_4 = sqrt( (2*N+1) / (2*N)  );

P(1:N-1) = 1./(N-m) .* ( (2*N-1).*sinphi.*Pn_1(1:N-1).*chi_1 - ( N + m -1).*Pn_2.*chi_2 );
P(N) =  sinphi *  Pn_1(N) *chi_3;
P(N+1) = cosphi * Pn_1(N) *chi_4;

Pn_2 = Pn_1;
Pn_1 = P;

end