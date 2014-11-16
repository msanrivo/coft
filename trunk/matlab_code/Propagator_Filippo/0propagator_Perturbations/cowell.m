%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% COWELL.M
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script propagates the state of an object in Earth orbit from 
% initial position and velocity vectors in cartesian coordinates using 
% Cowell method (see Vallado).  
%
% Author: Paloma Belen Mangas Carbajo      
%
% 
% Required input - in script -:
%
%       t0      [1x1]   initial time
%       tf      [1x1]   final time
%	x0      [6x1]   initial state. Initial  position and velocity 
%                       vectors in cartesian coordinates
%       pforces         vector of function handles to perturbation forces 
%       		   
% Output:
%
%       xOut    [6xn]   State from t0 to tf in n instants
% 
%
% Version: 1.0 28/06/2013
%          2.0 10/08/2013 Remove screen inputs
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters. Integration
%
tol     = 1e-12;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
%
% Input
% 
x0       = zeros([6,1]);
x0(1)    = 7500; 
x0(5)    = 7.2902*cos(30); 
x0(6)    = 7.2902*sin(30); 
t0       = 0; 
tf       = 86400; 
tspan    = [t0,tf];
%
pforces  = {@J2}; 
%
% Call to integrator
%
[t,xOut] = ode45( @perturb , tspan , x0 , options , pforces );
%
%
% Call to posprocess
%
