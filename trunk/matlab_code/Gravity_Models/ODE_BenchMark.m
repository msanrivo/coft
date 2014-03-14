clc
clear all
close all

%{
    Integrators BenchMark

%}

load('itggoce02_cell.mat');
[adaasdf, NN] = size(itggoce02_cell);

%% Declaring Physical Constants

GM=3986004.415e+8; % m^3/s^2
R = 6378136.3; %m

%% Declaring Initial Position

r = 1.25*R;
phi = 0.25;
lambda = 0;
[x,y,z] = sph2cart(lambda,phi,r);

%% Initial Velocity

vy    = 7500;
vx    = 7.2902*cos(30);
vz    = 7.2902*sin(30);

%% Calling Integrator
tf = 10;

tspan = [0 tf];
init_vals = [x;y;z;vx;vy;vz];

elements = 30;

sdom = floor(linspace(2,NN,elements));
time_matrix = nan(elements,5);

tol     = 1e-12;
opt = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);

for ct = 1:elements


    fprintf('Progress %1.3f %',ct/elements*100);
    
    tic
    [t,X] = ode113(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,sdom(ct),tf);
    time_matrix(ct,1) = toc;

    tic
    [t,X] = ode15s(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,sdom(ct),tf);
    time_matrix(ct,2) = toc;

    tic
    [t,X] = ode45(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,sdom(ct),tf);
    time_matrix(ct,3) = toc;

    tic
    [t,X] = ode23(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,sdom(ct),tf);
    time_matrix(ct,3) = toc;

    tic
    [t,X] = ode23s(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,sdom(ct),tf);
    time_matrix(ct,4) = toc;

    tic
    [t,X] = ode23t(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,sdom(ct),tf);
    time_matrix(ct,4) = toc;

end

plot(time_matrix(:,1),sdom)
hold on
plot(time_matrix(:,2),sdom,'k')
hold on
plot(time_matrix(:,3),sdom,'r')
hold on
plot(time_matrix(:,4),sdom,'g')
hold on
plot(time_matrix(:,5),sdom,'m')
