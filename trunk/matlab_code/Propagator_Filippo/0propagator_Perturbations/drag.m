% Calcula el efecto perturbador de la resistencia atmosferica.
%
% function dd_dt = drag( x, r ) 
% 
% Autor: Paloma Belen Mangas Carbajo     
 
% Entrada:
%       
%	x, vector posición/velocidad
%	r, sqrt(x(1)^2+x(2)^2+x(3)^2)
%          
% Salida:
%
%       dd_dt, vector con la v_pert y la a_pert
%
%           dd_dt(1), v_pert(1)
%           dd_dt(2), v_pert(2)
%           dd_dt(3), v_pert(3)
%           dd_dt(4), a_pert(1)
%           dd_dt(5), a_pert(2)
%           dd_dt(6), a_pert(3)

function [ dd_dt ] = drag(x, r)

Cd=2;
A=3.6; %[m^2]
m=1350; %[kg] 
H=200000; %[m]
r0=7298145; %[m]

rho=4e-13; %[kg/m^3]
w=7.29211585530066e-5; %[rad/s]
vr=sqrt((x(4)+w*x(2))^2+(x(5)-w*x(1))^2+x(6)^2);
rha=rho*exp((-r+r0)/H);

dd_dt(1)=0;
dd_dt(2)=0;
dd_dt(3)=0;
dd_dt(4)=-(1/2)*Cd*(A/m)*rha*vr*(x(4)+w*x(2));
dd_dt(5)=-(1/2)*Cd*(A/m)*rha*vr*(x(5)-w*x(1));
dd_dt(6)=-(1/2)*Cd*(A/m)*rha*vr*x(6);

return