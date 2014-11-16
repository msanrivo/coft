% function um_tf_xo_cowell 
% 
% Autor:   Paloma Belen Mangas Carbajo
 
% Entrada:
%
%   f, TLE, COE, ECI  
%   tf, tiempo de simulacion
%   x0, vector posición/velocidad
%   perturb_J2, activa/desactiva esta perturbacion
%   perturb_drag, activa/desactiva esta perturbacion
%   perturb_flux, activa/desactiva esta perturbacion
%	   
% Salida:
%
%   q, tiempo de simulacion / um 
%   x0, vector posición/velocidad
%   perturb_J2, activa/desactiva esta perturbacion
%   perturb_drag, activa/desactiva esta perturbacion
%   perturb_flux, activa/desactiva esta perturbacion

function [x0,q,perturb_J2,perturb_drag,perturb_flux] = um_tf_xo_cowell

G=6.67e-11; %en [m^3/s^2*kg]
M=5.9760179910044977511244377811094e24; %en [kg]

format long g

[f,tf,x0,perturb_J2,perturb_drag,perturb_flux]=prueba;

if f==1;

line_1= input( 'line_1: EXAMPLE [1 25544 98067    08264.51782528  -.00002182  00000-0 -11606-4 0  2927] ');
line_2= input( 'line_2: EXAMPLE [2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537] ');

S=1;
N=0;

g=input ( '¿Desea calcular h/d e incremento de t? [S/N] ');

if g==1;

tf=dia_hora(line_1);

elseif g==0;

tf=input( 'Introduzca tiempo total de simulacion: en [s] ');

end

q(1)=G;
q(2)=M;
q(3)=tf; 
b=line_2(8)*2*pi/86400;
um=q(1)*q(2);
a_3=um/b^2;
a=nthroot(a_3,3);
p0=[a 1e-7*line_2(5) line_2(3) line_2(4) line_2(6) line_2(7)]; %en [m, grados]

[r,v] = coe2rv(p0);

x0=[r(1) r(2) r(3) v(1) v(2) v(3)];

elseif f==2;

q(1)=G;
q(2)=M;
q(3)=tf; 

elseif f==3;

p0=input( 'Introduzca [a e i omega w M]: en [m, grados]' ); 

[r,v] = coe2rv(p0);

x0=[r(1) r(2) r(3) v(1) v(2) v(3)];

tf=input( 'Introduzca tiempo total de simulacion: en [s] ');

q(1)=G;
q(2)=M;
q(3)=tf; 

end