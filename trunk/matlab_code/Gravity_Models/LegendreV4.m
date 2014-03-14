function [P, Pn_1, Pn_2]  = LegendreV4(N,Pn_1,Pn_2,sinphi,cosphi)

%% n >= 2

P = nan(N+1,1);
m = [0:N-2]';

chi_1 = sqrt(  (2*N+1).*(N-m) ./( (2*N-1).*(N+m) )  );
chi_2 = sqrt(  (2*N+1).*(N-m).*(N-m-1) ./ ( (2*N-3).*(N+m).*(N+m-1) )  );

mm = N-1;
chi_3 = sqrt(2*mm+1);
chi_4 = sqrt( (2*N+1) / (2*N)  );

P(1:N-1) = 1./(N-m) .* ( (2*N-1).*sinphi.*Pn_1(1:N-1).*chi_1 - ( N + m -1).*Pn_2.*chi_2 );
P(N) =  sinphi *  Pn_1(N) *chi_3;
P(N+1) = cosphi * Pn_1(N) *chi_4;

Pn_2 = Pn_1;
Pn_1 = P;

end