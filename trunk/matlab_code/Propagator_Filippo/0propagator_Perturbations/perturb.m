%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% PERTURB.M
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function dxdt = perturb( t, x, pforce ) 
% 
% This function computes the times derivatives of the state including 
% the main gravitational term and acceleration due to perturbation forces
% 
%
% Author: Paloma Belen Mangas Carbajo      
%
% 
% Input:
%
%       t       [1x1]   time, independent variable
%	x       [6x1]   state: position and velocity vector in cartesian
%       pforce          vector of handles of the perturbative forces
%       		   
% Output:
%
%       pacc    [3x1]   perturbing acceleration vector in cartesian
% 
%
% Version: 1.0 28/06/2013
%          2.0 10/08/2013 Change in input and output of the functions
%                         vector of handles instead of flags 
%
%
function dxdt = perturb( t, x, pforce )
%
% Prameters
%
G       = 6.67e-20;  % [km^3/s^2*kg]
Mearth  = 5.9760179910044977511244377811094e24; % [kg]
um      = G*Mearth; 
%
% Initialization.
%
pacc = zeros([3,1]); 
dxdt = zeros([6,1]); 
%
% Loop: call to functions that compute perturbing forces 
%
for ii = 1:length( pforce )
	   pacc = pacc + pforce{ii}( t , x );
end
%
% Time derivatives
%
dxdt(1:3) = x(4:6); 
dxdt(4:6) = - um * x(1:3) / norm(x(1:3))^3 + pacc; 
%
return
