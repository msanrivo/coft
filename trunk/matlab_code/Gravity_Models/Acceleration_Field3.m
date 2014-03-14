function rdotdot = Acceleration_Field3(r,phi,lambda,GM,R,itggoce02_cell,NNN)

%{
    rdotdot is the output vector containing the acceleration at x,y,z in
    the vector form [x;y;z]

    Two degree counters are defined suchthat n is the degree for V and W
    whereas N is the degree of d2r_{N,m} (acceleration vector)

%}

%% Required Data

sinphi = sin(phi);
cosphi = cos(phi);

d2x = 0;
d2y = 0;
d2z = 0;

%% Computing Accelerations at n = 0

n = 0;
m = 0;

% Required Cnm @ n

Cnm = itggoce02_cell{n+1}(:,1);
Snm = itggoce02_cell{n+1}(:,2);

% Required Functions V & W @(n+1)

P = 1; %Legendre @ n = 0
Pn_1 = P;

P(1,1) = sinphi * Pn_1 * sqrt(3); %Legendre @ n = 1
P(2,1) = cosphi * Pn_1 * sqrt(3); %Legendre @ n = 1

Pn_2 = Pn_1;
Pn_1 = P;

M = [0;1];

V = (R/r)^(n+1) .* cos(M.*lambda) .* P; % Normalized
W = (R/r)^(n+1) .* sin(M.*lambda) .* P; % Normalized

NF = sqrt(1/3);

d2x = d2x - Cnm(1)*V(2) * NF;
d2y = d2y - Cnm(1)*W(2) * NF;
d2z = d2z + (-Cnm(1)*V(1) - Snm(1)*W(1)) * NF;

%% Computing Accelerations at n = 1

n = n+1; %n = 1

% Required Cnm @ n

Cnm = itggoce02_cell{n+1}(:,1);
Snm = itggoce02_cell{n+1}(:,2);

% Required Functions V & W @(n+1)

[P, Pn_1, Pn_2]  = LegendreV4(n+1,Pn_1,Pn_2,sinphi,cosphi);
M = [0;1;2];
V = (R/r)^(n+1) .* cos(M.*lambda) .* P; % Normalized
W = (R/r)^(n+1) .* sin(M.*lambda) .* P; % Normalized

%%% m = 0;

NF_x = sqrt(0.5*(2*n+1)/(2*n+3)*(n+2)*(n+1));
NF_y = NF_x;
NF_z = sqrt((2*n+1)/(2*n+3)*(n+1));

d2x = d2x - Cnm(1)*V(2) * NF_x;
d2y = d2y - Cnm(1)*W(2) * NF_y;
d2z = d2z + (n+1)*(-Cnm(1)*V(1) - Snm(1)*W(1)) * NF_z;

%%% m = 1;
m = 1;

NF_1 = sqrt((2*n+1)/(2*n+3)*(n+m+2)*(n+m+1));
NF_2 = sqrt(2*(2*n+1)/(2*n+3)/(n-m+2)/(n-m+1)); % Notice sqrt(2) is an exception just at n = 1
NF_3 = sqrt((2*n+1)/(2*n+3)*(n+m+1)/(n-m+1));

d2x = d2x + 0.5*((-Cnm(2)*V(3) - Snm(2)*W(3))* NF_1 + (n-m+2)*(n-m+1)*(Cnm(2)*V(1) + Snm(2)*W(1))*NF_2 );
d2y = d2y + 0.5*((-Cnm(2)*W(3) + Snm(2)*V(3))* NF_1 + (n-m+2)*(n-m+1)*(-Cnm(2)*W(1) + Snm(2)*V(1))*NF_2 );
d2z = d2z + (n-m+1)*(-Cnm(2)*V(2) - Snm(2)*W(2)) * NF_3;


%% Looping the function from n = 2 : NN

for n = 2 : NNN-1

    %% Data Required @n

    Cnm = itggoce02_cell{n+1}(:,1);
    Snm = itggoce02_cell{n+1}(:,2);

    % Required Functions V & W @(n+1)

    [P, Pn_1, Pn_2]  = LegendreV4(n+1,Pn_1,Pn_2,sinphi,cosphi);
    M = [0:(n+1)]';
    V = (R/r)^(n+1) .* cos(M.*lambda) .* P; % Normalized
    W = (R/r)^(n+1) .* sin(M.*lambda) .* P; % Normalized

    % Computing @ m = 0

    NF_x = sqrt(0.5*(2*n+1)/(2*n+3)*(n+2)*(n+1));
    NF_y = NF_x;
    NF_z = sqrt((2*n+1)/(2*n+3)*(n+1));

    d2x = d2x - Cnm(1)*V(2) * NF_x;
    d2y = d2y - Cnm(1)*W(2) * NF_y;
    d2z = d2z + (n+1)*(-Cnm(1)*V(1) - Snm(1)*W(1)) * NF_z;

    % Computing @ m = 1
    m = 1;

    NF_1 = sqrt((2*n+1)/(2*n+3)*(n+m+2)*(n+m+1));
    NF_2 = sqrt(2*(2*n+1)/(2*n+3)/(n-m+2)/(n-m+1)); % Notice sqrt(2) is an exception just at n = 1
    NF_3 = sqrt((2*n+1)/(2*n+3)*(n+m+1)/(n-m+1));

    d2x = d2x + 0.5*((-Cnm(2)*V(3) - Snm(2)*W(3))* NF_1 + (n-m+2)*(n-m+1)*(Cnm(2)*V(1) + Snm(2)*W(1))*NF_2 );
    d2y = d2y + 0.5*((-Cnm(2)*W(3) + Snm(2)*V(3))* NF_1 + (n-m+2)*(n-m+1)*(-Cnm(2)*W(1) + Snm(2)*V(1))*NF_2 );
    d2z = d2z + (n-m+1)*(-Cnm(2)*V(2) - Snm(2)*W(2)) * NF_3;
    
    % Computing @ m > 1

    for m = 2:n

        NF_1 = sqrt((2*n+1)/(2*n+3)*(n+m+2)*(n+m+1));
        NF_2 = sqrt((2*n+1)/(2*n+3)/(n-m+2)/(n-m+1)); % Notice sqrt(2) is an exception just at n = 1
        NF_3 = sqrt((2*n+1)/(2*n+3)*(n+m+1)/(n-m+1));

        d2x = d2x + 0.5*((-Cnm(m+1)*V(m+2) - Snm(m+1)*W(m+2))* NF_1 + (n-m+2)*(n-m+1)*(Cnm(m+1)*V(m) + Snm(m+1)*W(m))*NF_2 );
        d2y = d2y + 0.5*((-Cnm(m+1)*W(m+2) + Snm(m+1)*V(m+2))* NF_1 + (n-m+2)*(n-m+1)*(-Cnm(m+1)*W(m) + Snm(m+1)*V(m))*NF_2 );
        d2z = d2z + (n-m+1)*(-Cnm(m+1)*V(m+1) - Snm(m+1)*W(m+1)) * NF_3;

    end

end

rdotdot = GM/R^2 * [d2x;d2y;d2z];

end