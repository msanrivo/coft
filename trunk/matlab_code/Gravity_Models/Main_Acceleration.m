clc
clear all
close all
format long

load('itggoce02_cell.mat');
[adaasdf, NN] = size(itggoce02_cell);

%% Declaring Physical Constants

GM=3986004.415e+8; % m^3/s^2
R = 6378136.3; %m

%% Declaring Initial Position

r = 1.35*R;
phi = pi/2;
lambda = 0;
[x,y,z] = sph2cart(lambda,phi,r);

%% Initial Velocity

vy    = 7500; 
vx    = 7.2902*cos(30); 
vz    = 7.2902*sin(30); 

%% Calling Integrator
tf = 86400; % to see deviation due to spherical harmonics... huge time required

tspan = [0 tf];
init_vals = [x;y;z;vx;vy;vz];

tol     = 1e-12;
opt = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);

tic
[t,X] = ode113(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,NN,tf);
toc


%% In Case of a long computation may be interesting to save results

%%%%% I have saved onced for tspan [0 86400*10] seconds at tol = 1e-12;

% save('time_vector','t');
% save('State_Vector','X');


%% Simulation Running

[xx,yy,zz] = sphere(100);
xx = xx.*R;
yy = yy.*R;
zz = zz.*R;

for i = 1 : length(t)
   
    surf(xx,yy,zz);
    hold on
    plot3(X(:,1),X(:,2),X(:,3));
    hold on
    plot3(X(i,1),X(i,2),X(i,3),'r+');
    hold off
    title(sprintf('Time passed: %1.1f', t(i)));
    
    drawnow;
    
    
end

