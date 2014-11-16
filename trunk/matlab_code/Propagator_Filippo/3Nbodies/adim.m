%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script adimensionalizes the initial conditions and time span
%
% Author: Pablo A. Machuca Varela
%
% Input:
% b          [1x1] Chosen body to determine characteristic magnitudes
% M          [nx1] Masses of the n bodies
% mu         [nx1] Gravitational parameters
% ro         [1x3*n]    [xo, yo, zo]      Initial positions
% vo         [1x3*n]    [vxo, vyo, wvzo]  Initial velocities
% tspan
%
% Output:
% mc         [1x1]  Characteristic mass
% Lc         [1x1]  Characteristic length
% tc         [1x1]  Characteristic time
% MA         [nx1]  Adimensionalized masses of the n bodies
% muA        [nx1]  Adimensionalized gravitational parameters
% roA        [1x3*n]    [xoA, yoA, zoA]      Adimensionalized initial positions
% voA        [1x3*n]    [vxoA, vyoA, vzoA]   Adimensionalized initial velocities
% tspanA     Adimensionalized tspan
% InitCond   [1x6*n]'   [xoA, yoA, zoA, vxoA, vyoA, vzoA]' Adimensionalized initial positions and velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Chosen body
b=7;
%
%Characteristic magnitudes
mc = M(b); 
Lc = norm(ro((3*(b-1)+1):(3*(b-1)+3))); %norm([x, y, z])
tc = 2*pi*norm(ro((3*(b-1)+1):(3*(b-1)+3)))^(3/2)/((mu(1))^(1/2)); %Third Kepler's Law: T^2 = 4*pi^2/(G*M)*r^3
%
%Adimensionalization
MA  = M/mc;
muA = mu/(Lc^3/tc^2);
roA = ro/Lc;
voA = vo/(Lc/tc);
tspanA = tspan/tc;
%
InitCond = [roA,voA]';
