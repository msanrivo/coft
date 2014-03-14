clc
clear all
close all

%{

    Testing GeoPotential Convergence at a single point in space arbitrary chosen
    and with the given values of physical constants

%}

%% Declaring Physical Constants

GM=3986004.415e+8; % m^3/s^2
R = 6378136.3; %m

%% Loading Coefficients DataBase

load('itggoce02_cell.mat');
[aa, NN] = size(itggoce02_cell);

%% Declaring Position

for i = 1:5

r = (1+rand)*R;
phi = rand; sinphi = sin(phi); cosphi = cos(phi);
lambda = rand;

dim = 200;
dom = floor(linspace(1,NN-1,dim));
U = nan(dim,1);
ct = 1;


for ss = dom
    
    outsum = 0;

    %% Legendre @ n = 1

    Pn_1 = 1;
    Pn_2 = 0;

    P(1,1) = sinphi * Pn_1 * sqrt(3);
    P(2,1) = cosphi * Pn_1 * sqrt(3);

    Cnm = itggoce02_cell{2}(:,1);
    Snm = itggoce02_cell{2}(:,2);
    m = [0:1]';
    outsum = outsum + sum( (R/r)^1 .* P .* ( Cnm.*cos(m.*lambda) + Snm.*sin(m.*lambda) ) );

    Pn_2 = Pn_1;
    Pn_1 = P;

    for n = 2:ss

        Cnm = itggoce02_cell{n+1}(:,1);
        Snm = itggoce02_cell{n+1}(:,2);

        m = [0:n]';

        [P, Pn_1, Pn_2]  = LegendreV4(n,Pn_1,Pn_2,sinphi,cosphi);

        outsum = outsum + sum( (R/r)^n .* P .* ( Cnm.*cos(m.*lambda) + Snm.*sin(m.*lambda) ) );

    end
    
    clear P
    clear Pn_1
    clear Pn_2
    
    U(ct) = (1+outsum)*GM/r;
    ct = ct+1;

end

semilogy(dom(1:end-1),abs(diff(U)),'Color',rand(1,3));
hold on
axis([0 NN 0 max(abs(diff(U)))*10]);
xlabel('Model Size')
ylabel('GeoPotential Difference \SigmaU_{0}^{n} -\SigmaU_{0}^{n-1}');
title('GeoPotential Convergence Test');

% legend(sprintf('( r/R, \\phi, \\lambda) = (%1.3f, %1.3f, %1.3f)',r/R,phi*180/pi,lambda*180/pi));
end