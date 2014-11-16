function U = GeoPotential_Calculator(GM,R,r,phi,lambda,Geop_Coeff,N)

%{

    GeoPotential Computation Function Based On Spherical Harmonics

    Author: Diego García Pardo (UNIVERSIDAD CARLOS III DE MADRID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Inputs:
        
        -- GM:          Standard Gravitational Parameter        [m^3/s^2]
        -- R            Planet Radious                          [m]
        -- r:           Distance to Earth CM (r>R)              [m]
        -- phi:         Latitude                                [rad]
        -- lambda:      longitude                               [rad]
        -- Geop_Coeff:  Spherical Harmonics Coeffs              []
        -- N:           Model Size Used                         []

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Output:

        -- U  =          Geopotential Value at the given location [J]
        -- dU =          Geopotential Value at the given location [J]


%}


%% Seed Legendre Polynomial

sinphi = sin(phi);
cosphi = cos(phi);
U = 1;

Pn_1 = 1;                                                   % Legendre Polynomial @ n = 0
P = [sinphi * Pn_1 * sqrt(3); cosphi * Pn_1 * sqrt(3)];     % Legendre Polynomial @ n = 1

Cnm = Geop_Coeff{2}(:,1);
Snm = Geop_Coeff{2}(:,2);

m = (0:1)';
U = U + sum( (R/r)^1 .* P .* ( Cnm.*cos(m.*lambda) + Snm.*sin(m.*lambda) ) );
Pn_2 = Pn_1;
Pn_1 = P;

if N >= 2
    for n = 2:N
        
        Cnm = Geop_Coeff{n+1}(:,1);
        Snm = Geop_Coeff{n+1}(:,2);
        
        [P, Pn_1, Pn_2]  = Normalized_Legendre(n,Pn_1,Pn_2,sinphi,cosphi);
        m = (0:n)';
        
        U = U + sum( (R/r)^n .* P .* ( Cnm.*cos(m.*lambda) + Snm.*sin(m.*lambda) ) );
        
    end
end

U = GM/r * U;

end
