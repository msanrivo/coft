%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script dimensionalizes the results
%
% Author: Pablo A. Machuca Varela
%
% Input:
% T          [length(T)x1]    Adimensionalized time
% X          [length(T)x6*n]  [x, y, z, vx, vy, vz] Adimensionalized results from the integration
% MA         [nx1]  Adimensionalized masses of the n bodies
% muA        [nx1]  Adimensionalized gravitational parameters
% mc         [1x1]  Characteristic mass
% Lc         [1x1]  Characteristic length
% tc         [1x1]  Characteristic time
%
% Output:
% mcR        [1x1]   Characteristic mass in the output units
% LcR        [1x1]   Characteristic length in the output units
% tcR        [1x1]   Characteristic time in the output units
% M          [nx1]   Dimensionalized masses of the n bodies
% mu         [nx1]   Dimensionalized gravitational parameters
% T          [length(T)x1]    Dimensionalized time
% X          [length(T)x6*n]  [x, y, z, vx, vy, vz] Dimensionalized results from the integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%1 AU=149597870691m 1d=86400s
%
%Characteristic magnitudes multiplied by the conversion factors
mcR = mc*1; 
LcR = Lc*149597870691; 
tcR = tc*86400; 

M  = MA*mcR;
mu = muA*(LcR^3/tcR^2);
X(:,1:3*n) = X(:,1:3*n)*LcR;
X(:,(3*n+1):end) = X(:,(3*n+1):end)*(LcR/tcR);
T = T*tcR;
