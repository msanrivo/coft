clc
clear all
close all

%{

    ORBIT PROPAGATOR BASED ON SPHERICAL HARMONICS

    AUTHOR: Diego García Pardo (UNIVERSIDAD CARLOS III DE MADRID)

    This Script Makes Use of the functions:
    
        --  Acc_Field In order to obtain the inertial Acceleration

        --  FictAcc in order to obtain the Fictitous forces related to a
            frame spinning with Earth and whose origin is at Earth CM
        
        --  Field_Integrator In order to merge kinematic and dynamic
            relationships

%}

%% Defining Planet Characteristics

GM= 3.9860044150e+14;                   % m^3/s^2 
R = 6.3781363000e+06;                   % m

load('itggoce02_cell_modified.mat');
[~, NN] = size(itggoce02_cell);         % Model Size Evaluation


%% Declaring Initial Conditions For Integrator

%%% Initial Position

r = 1.1*R;                              % Distance From Earth CM
phi = pi/4;                             % Planet Latitude           (Geocentric)
lambda = pi/5;                             % Planet Longitude          (Geocentric)
[x,y,z] = sph2cart(lambda,phi,r);       % Cartesian Transformation

%%% Initial Velocity

vx    = 7500;                           % [m/s]
vy    = 7.2902*cos(30);                 % [m/s]
vz    = 7.2902*sin(30);                 % [m/s]

%% Integrator Set-Up

t0 = 0; tf = 1e-4; 
tspan = [t0 tf];                        % Interval Of Time Of Integration [secs]

tol     = 1e-12;
opt = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);

%%% Integrator Call (NOTE: INTEGRATION IS DONE IN CARTESIAN COORDINATES)

tic

[t,X] = ode113(@ Field_Integrator,tspan,[x;y;z;vx;vy;vz],opt,GM,R,itggoce02_cell,NN,tf);

T = toc;                                % Computational Time Spend

fprintf('\nIt took %1.3f [secs] to Compute %1.3f Real Time Simulation',T,tf-t0)
fprintf('\nIt takes %1.3f [secs] to compute one real second (Simulation Average)\n',T/(tf-t0));


%% Results Representation

Simulation(t,X,R);
