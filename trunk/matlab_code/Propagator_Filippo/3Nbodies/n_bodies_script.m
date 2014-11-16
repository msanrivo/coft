%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script defines the integration time and calls the other functions to run the full code
%
% Author: Pablo A. Machuca Varela
%
% Input:
% octave    [1x1]   Possible values: 0, 1. To define whether the code is run in Matlab: 0, or in Octave: 1
% tspan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
close all
clc
%%%%%%%%%%%%
octave = 0;
%%%%%%%%%%%%
%
n_bodies_preprocess %It calls the function where initial conditions are set
%
tspan = [0 (2*365)];
%
adim %It calls the function where initial conditions are adimensionalized
%%%%%%%%%%%%%%
%Integration
%
% Octave
if (octave)
tic
tspan_Oc = linspace(tspanA(1), tspanA(end)); 
X = lsode(@(x,t)n_bodies_function_Oc(x,t,muA,n), InitCond, tspan_Oc);
toc
T = tspan_Oc; 
else
%
% Matlab
options=odeset('AbsTol',1e-14,'Reltol',1e-11);
tic
[T X]=ode45(@(t,x)n_bodies_function(t,x,muA,n), tspanA, InitCond, options);
toc
end
%%%%%%%%%%%%%%
%
redim %It calls the function where data is dimensionalized
%
n_bodies_postprocess %It calls the function where the results are postprocessed
