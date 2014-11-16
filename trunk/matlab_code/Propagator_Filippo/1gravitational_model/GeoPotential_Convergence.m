clc
clear all
close all

%{

    Geopotential Convergence Analysis

    Author: Diego García Pardo (UNIVERSIDAD CARLOS III DE MADRID)


%}

%% Defining Planet Characteristics

GM= 3.9860044150e+14;                   % m^3/s^2
R = 6.3781363000e+06;                   % m

load('itggoce02_cell_modified.mat');
[~, N] = size(itggoce02_cell);          % Model Size Evaluation

%% Defining Location

rR = 1.1;  r = rR*R;
phidom = linspace(-pi/2,pi/2,100);
ldom = linspace(0,2*pi,100);

S = nan(length(phidom),length(ldom));
C = S;
P = S;
L = S;
UU = S;

for i = 1:length(phidom)
    
    for j = 1:length(ldom)
        
        display([i,j])
        
        P(i,j) = phidom(i);
        L(i,j) = ldom(j);
        
        
        U = nan(N,1);
        U(1) = GeoPotential_Calculator(GM,R,r,phidom(i),ldom(j),itggoce02_cell,1);
        
        for n = 2:N
            
            U(n) = GeoPotential_Calculator(GM,R,r,phidom(i),ldom(j),itggoce02_cell,n);
            
            if abs(U(n) - U(n-1))  < 1e-24
                
                S(i,j) = n;                         % Model Size Used
                C(i,j) = min(abs(diff(U(1:n-1))));       % Convergence Acquired
                UU(i,j) = U(n);
                break
                
            end
            
            %             semilogy(abs(diff(U)));
            %             title(sprintf('Location [r/R, \\phi, \\lambda] = [%1.1f, %1.3f º, %1.3f º]',rR,P(i,j),L(i,j)));
            %             drawnow
            
        end
        
        S(i,j) = n;
        C(i,j) = min(abs(diff(U(1:n-1))));
        UU(i,j) = U(n);
        
    end
    
end

%% World Convergence And Model Size
load coast

for i = 1:length(phidom)
    for j =1:length(ldom)
    
        P(i,j) = phidom(i);
        L(i,j) = ldom(j);
        
    end
end


P = P*180/pi;
L = L*180/pi -180; %-180 º to correct location in map (function geoshow)
UU = UU./1e6;

%%
close all

% DATA{1} = P;
% DATA{2} = L;
% DATA{3} = UU;
%  
% save('Geopotential_Convergence_Results.mat','DATA');


figure
subplot(2,1,1)
contour(L,P,log10(C)); colorbar;
geoshow(lat,long,'Color','black','LineWidth',3)
title(sprintf('EQ: 1.25, Relative Convergence At constant Radious r/R = %1.3f',rR));
xlabel('Longitude [deg] ');
ylabel('Latitude [deg] ');

subplot(2,1,2)
contour(L,P,(S)); colorbar;
geoshow(lat,long,'Color','black','LineWidth',3)
title(sprintf('Model Size Employed To Reach Convergence At constant Radious r/R = %1.3f',rR));
xlabel('Longitude [deg] ');
ylabel('Latitude [deg] ');

figure
contourf(L,P,UU); colorbar;
geoshow(lat,long,'Color','Black')
title(sprintf('EQ: 1.26, Potential Energy Levels At constant Radious r/R = %1.3f \n Contours in MegaJoules',rR));
xlabel('Longitude [deg] ');
ylabel('Latitude [deg] ');

