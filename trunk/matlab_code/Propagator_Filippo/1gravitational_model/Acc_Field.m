function rdotdot = Acc_Field(r,phi,lambda,GM,R,GeoP_Coeff,NN)
%{

    Computation Of Inertial Acceleration Field Associtaed with Planetary Characteristics:

    Author: Diego García Pardo (UNIVERSITY CARLOS III OF MADRID)

        --Std. Gravitational Constant       GM  [m^3 s^-2]
        --Planet Radius                     R   [m]
        --Space Location
                - Distance from Planet CM   r   [m]
                - Longitude                 phi [rads]
                - Latitude                  lambda [rads]

        --Geopotential Coefficients (Fully Normalized)
                - Cell Size of NNx2, 1st column Cnm 2nd Snm
                
    It is used the scheme shown in Montebruck (Satellite Orbits)

    NOTE: This Algorithm requires NN >= 2

%}

%% Required Parameters in Algorithm

sinphi = sin(phi);
cosphi = cos(phi);


%% Order 0

n = 0;              % Model Order Counter;
m = 0;              % Associated Index

N = n+1;            % Counter for functions V & W (which are 1 order higher)
M = (0:1:N)';      % Associated Index for functions V & W

%%% Zonal Coefficients Associated Acceleration @(m = 0)

Cnm = GeoP_Coeff{n+1}(:,1);
Snm = GeoP_Coeff{n+1}(:,2);

P00 = 1;
P1m = [sinphi*sqrt(3) ; cosphi*sqrt(3)];

V = (R/r)^(N+1) * P1m .* cos(M.*lambda);
W = (R/r)^(N+1) * P1m .* sin(M.*lambda);

NF1  = sqrt(0.5*(2*n+1)*(n+2)*(n+1)/(2*n+3));
NF2  = sqrt((2*n+1)*(n+m+1)/(2*n+3)/(n-m+1));

d2x = -Cnm * V(2) * NF1;
d2y = -Cnm * W(2) * NF1;
d2z = (n+1)*(-Cnm*V(1) - Snm*W(1)) * NF2;

%% Counter increase to order 1

n = n+1;            % n = 1;
% m = 0:1:n;

N = n+1;
M = (0:1:N)';

Cnm = GeoP_Coeff{n+1}(:,1);
Snm = GeoP_Coeff{n+1}(:,2);

[P,Pn_1,Pn_2] = Normalized_Legendre(N,P1m,P00,sinphi,cosphi);

V = (R/r)^(N+1) * P .* cos(M.*lambda);
W = (R/r)^(N+1) * P .* sin(M.*lambda);

%%% Zonal Coefficients Associated Acceleration @(m = 0)
m = 0;
NF1  = sqrt(0.5*(2*n+1)*(n+2)*(n+1)/(2*n+3));
NF2  = sqrt((2*n+1)*(n+m+1)/(2*n+3)/(n-m+1));

d2x = -Cnm(1) * V(2) * NF1                              + d2x;
d2y = -Cnm(1) * W(2) * NF1                              + d2y;
d2z = (n-m+1)*(-Cnm(1)*V(1) - Snm(1)*W(1)) * NF2        + d2z;


%%% Acceleration @(m = 1), Scheme Exception Due To Normalization
m = 1;
NF1 = sqrt((2*n+1)/(2*n+3)*(n+m+2)*(n+m+1));
NF2 = sqrt(2*(2*n+1)/(2*n+3)/(n-m+2)/(n-m+1));
NF3 = sqrt((2*n+1)*(n+m+1)/(2*n+3)/(n-m+1));

d2x = d2x   +   0.5*( (-Cnm(2)*V(3) - Snm(2)*W(3)) * NF1 + (n-m+2)*(n-m+1)*(+Cnm(2)*V(1) + Snm(2)*W(1)) * NF2 ) ;
d2y = d2y   +   0.5*( (-Cnm(2)*W(3) + Snm(2)*V(3)) * NF1 + (n-m+2)*(n-m+1)*(-Cnm(2)*V(1) + Snm(2)*V(1)) * NF2 ) ;
d2z = d2z   +   (n-m+1)*(-Cnm(2)*V(2) - Snm(2)*W(2)) * NF3;


for n = 2 : 1 : NN-1    % The algorithm only gets to NN-1 because it is a forward recursion so at NN-1 we compute the value at NN
    
    N = n+1;
    M = (0:1:N)';
    
    Cnm = GeoP_Coeff{n+1}(:,1);
    Snm = GeoP_Coeff{n+1}(:,2);
    
    [P,Pn_1,Pn_2] = Normalized_Legendre(N,Pn_1,Pn_2,sinphi,cosphi);
    
    V = (R/r)^(N+1) * P .* cos(M.*lambda);
    W = (R/r)^(N+1) * P .* sin(M.*lambda);
    
    %%% Zonal Coefficients Associated Acceleration @(m = 0)
    m = 0;
    NF1  = sqrt(0.5*(2*n+1)*(n+2)*(n+1)/(2*n+3));
    NF2  = sqrt((2*n+1)*(n+m+1)/(2*n+3)/(n-m+1));
    
    d2x = -Cnm(1) * V(2) * NF1                          + d2x;
    d2y = -Cnm(1) * W(2) * NF1                          + d2y;
    d2z = (n-m+1)*(-Cnm(1)*V(1) - Snm(1)*W(1)) * NF2    + d2z;
    
    %%% Acceleration @(m = 1), Scheme Exception Due To Normalization
    m = 1;
    NF1 = sqrt((2*n+1)/(2*n+3)*(n+m+2)*(n+m+1));
    NF2 = sqrt(2*(2*n+1)/(2*n+3)/(n-m+2)/(n-m+1));
    NF3 = sqrt((2*n+1)*(n+m+1)/(2*n+3)/(n-m+1));
    
    d2x = d2x   +   0.5*( (-Cnm(2)*V(3) - Snm(2)*W(3)) * NF1 + (n-m+2)*(n-m+1)*(+Cnm(2)*V(1) + Snm(2)*W(1)) * NF2 ) ;
    d2y = d2y   +   0.5*( (-Cnm(2)*W(3) + Snm(2)*V(3)) * NF1 + (n-m+2)*(n-m+1)*(-Cnm(2)*V(1) + Snm(2)*V(1)) * NF2 ) ;
    d2z = d2z   +   (n-m+1)*(-Cnm(2)*V(2) - Snm(2)*W(2)) * NF3;
    
    
    %%% Generality Of Terms
    for m = 2:1:n
        
        NF1 = sqrt((2*n+1)/(2*n+3)*(n+m+2)*(n+m+1));
        NF2 = sqrt((2*n+1)/(2*n+3)/(n-m+2)/(n-m+1));
        NF3 = sqrt((2*n+1)*(n+m+1)/(2*n+3)/(n-m+1));
        
        d2x = d2x   +  0.5*( (-Cnm(m+1)*V(m+2) - Snm(m+1)*W(m+2)) * NF1 + (n-m+2)*(n-m+1)*(+Cnm(m+1)*V(m) + Snm(m+1)*W(m)) * NF2 ) ;
        d2y = d2y   +  0.5*( (-Cnm(m+1)*W(m+2) + Snm(m+1)*V(m+2)) * NF1 + (n-m+2)*(n-m+1)*(-Cnm(m+1)*V(m) + Snm(m+1)*V(m)) * NF2 ) ;
        d2z = d2z   +  (n-m+1)*(-Cnm(m+1)*V(m+1) - Snm(m+1)*W(m+1)) * NF3;
        
        
    end
    
end

rdotdot = (GM/R^2) * [d2x;d2y;d2z];

end