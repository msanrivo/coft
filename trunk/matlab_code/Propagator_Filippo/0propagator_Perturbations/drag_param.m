% Calcula la variacion que se produce en los parametros orbitales debido a la resistencia atmosferica.
%
% function ddp_dt = drag_param( x, um )                                                                                                                                             
% 
% Autor: Paloma Belen Mangas Carbajo      
 
% Entrada:
%       
%       x, parametros orbitales
%       um, parametro gravitacional

%	   
% Salida:
%
%       ddp_dt, vector con las variaciones de los parametros
%
%	ddp_dt(1)=da_dt
%	ddp_dt(2)=de_dt  
%	ddp_dt(3)=di_dt 
%	ddp_dt(4)=dO_dt	
%	ddp_dt(5)=dw_dt 
%	ddp_dt(6)=dM_dt 

function [ ddp_dt ]=drag_param( x, um )

Cd=2;
A=3.6; %[m^2]
m=1350; %[kg]
H=200000; %[m]
r0=8118521.49; %[m]

rho=4e-13; %[kg/m^3]
rha=rho*exp((-x(1)+r0)/H);

ddp_dt(1)=-(Cd*A/m)*rha*sqrt(um*H/2*pi*x(2))*((1+x(2))/(1-x(2)))^2;
ddp_dt(2)=(1-x(2))*(-(Cd*A/m)*rha*sqrt(um*H/2*pi*x(2))*((1+x(2))/(1-x(2)))^2)/x(1);
ddp_dt(3)=0;
ddp_dt(4)=0;
ddp_dt(5)=0;
ddp_dt(6)=0;

return