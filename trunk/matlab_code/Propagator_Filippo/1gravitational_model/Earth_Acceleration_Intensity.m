clc
clear all
close all

%{

    Gravity On Earth Surface

    Author: Diego García Pardo (UNIVERSIDAD CARLOS III DE MADRID)


%}

%% Defining Planet Characteristics

GM= 3.9860044150e+14;                   % m^3/s^2
R = 6.3781363000e+06;                   % m

load('itggoce02_cell_modified.mat');
[~, N] = size(itggoce02_cell);          % Model Size Evaluation

%% Defining Location

rR = 1;  r = rR*R;
phidom = linspace(-pi/2,pi/2,100);
ldom = linspace(0,2*pi,30);

%% Computations

P = nan(length(phidom),length(ldom));
L = P;
A = P;

for j = 1:length(ldom)
    
    parfor i = 1:length(phidom)
        
        display([i,j])
        
        P(i,j) = phidom(i);
        L(i,j) = ldom(j);
        
        A(i,j) = norm(Acc_Field(r,phidom(i),ldom(j),GM,R,itggoce02_cell,N));
        
        
    end
    
end


%% World Map
load coast
P = P*180/pi;
L = L*180/pi -180; %-180 º to correct location in map (function geoshow)

%% Plotting Results
close all

contourf(L,P,A); colorbar; hold on;
geoshow(lat,long,'Color','black','LineWidth',2)
xlabel('Longitude [deg] ');
ylabel('Latitude [deg] ');
title(sprintf('Mean Acceleration On Earth Surface: %1.6f m s^-2',mean(mean(A))));
