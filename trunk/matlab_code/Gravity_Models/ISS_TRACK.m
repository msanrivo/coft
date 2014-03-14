clc
clear all
close all
format long

load('itggoce02_cell.mat');
[adaasdf, NN] = size(itggoce02_cell);

%% Declaring Physical Constants

GM=3986004.415e+8; % m^3/s^2
R = 6378136.3; %m

%% Declaring ISS Initial Position

r = 1.06758400569144*R;
phi = -0.771086463531095;
lambda = 2.87089208660547;
[x,y,z] = sph2cart(lambda,phi,r);

dphi = -0.000558020546957076;
dlambda = 0.00128184737285361;
dr = 10.1388888888889;

[vx,vy,vz] = cart_vel(r,phi,lambda,dr, dphi, dlambda);


%% Calling Integrator
tf = 3928; % to see deviation due to spherical harmonics... huge time required

tspan = [0 tf];
init_vals = [x;y;z;vx;vy;vz];

tol     = 1e-12;
opt = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);

tic
[t,X] = ode113(@ Field_Integrator,tspan,init_vals,opt,GM,R,itggoce02_cell,NN,tf);
toc

[th,phi,r] = cart2sph(X(end,1),X(end,2),X(end,3));
Longitude = th*180/pi
Latitude = phi*180/pi
Altitude = (r-R)/1000

%% Simulation Running

[xx,yy,zz] = sphere(100);
xx = xx.*R;
yy = yy.*R;
zz = zz.*R;

for i = 1 : length(t)
   
    subplot(1,2,1)
    mesh(xx,yy,zz);
    hold on
    plot3(X(:,1),X(:,2),X(:,3));
    hold on
    plot3(X(i,1),X(i,2),X(i,3),'r+');
    hold off
    title(sprintf('Time passed: %1.1f', t(i)));
    
    subplot(1,2,2)
    r = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2);
    plot(t,r./R);
    hold on
    plot(t,ones(length(t),1),'r');
    drawnow
        
    
end

