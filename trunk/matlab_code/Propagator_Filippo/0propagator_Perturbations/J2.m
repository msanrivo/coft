%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% J2.M
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function j2acc = J2( t, x) 
% 
% This function computes the acceleration due to the secular contribution
% of the Earth flattening 
%
% Author: Paloma Belen Mangas Carbajo      
%
% 
% Input:
%
%       t       [1x1]   time, independent variable
%	x       [6x1]   state: position and velocity vector in cartesian
%       		   
% Output:
%
%       j2acc    [3x1]   perturbing acceleration vector in cartesian
% 
%
% Version: 1.0 28/06/2013
%          2.0 10/08/2013 Change in input to make the function compatible  
%                         with perturb.m
%
%
function j2acc = J2( t, x) 
%
% Parameters
%
G       = 6.67e-20;  % [km^3/s^2*kg]
Mearth  = 5.9760179910044977511244377811094e24; % [kg]
Rearth  = 6378.145; % Equatorial radius of Earth [km]
um      = G*Mearth; 
J2      = 0.00108248;
%
% Variables
%
r  = norm(x(1:3)); % modulus of the position vector
r3 = r*r*r; 
%
% Initialization
%
j2acc = zeros([3,1]); 
%
% Computation
%
j2acc(1) = -(um*x(1)/r3)*(-(3/2)*J2*(Rearth/r)^2*(5*(x(3)/r)^2-1));
j2acc(2) = x(2)/x(1)*j2acc(1);
j2acc(3) = -(um*x(3)/r3)*(-(3/2)*J2*(Rearth/r)^2*(5*(x(3)/r)^2-3));
%
return
