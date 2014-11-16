% Calcula el efecto perturbador del achatamiento terrestre.
%
% function dj_dt = J2( x, um, R, r3 ) 
% 
% Autor: Paloma Belen Mangas Carbajo    
   
 
% Entrada:
%       
%   x, vector posición/velocidad
%   um, parametro gravitacional
%   R, Radio de la Tierra
%   r, r=(norm(x(1:3)))
%   r3, sqrt(x(1)^2+x(2)^2+x(3)^2)^3
%          
% Salida:
%
%       dj_dt, vector con la v_pert y la a_pert
%
%           dj_dt(1), v_pert(1)
%           dj_dt(2), v_pert(2)
%           dj_dt(3), v_pert(3)
%           dj_dt(4), a_pert(1)
%           dj_dt(5), a_pert(2)
%           dj_dt(6), a_pert(3)

function [ dj_dt ] = J2( x, um, R, r, r3)

J2=0.00108248;

dj_dt(1)=0;
dj_dt(2)=0;
dj_dt(3)=0;
dj_dt(4)=-(um*x(1)/r3)*(-(3/2)*J2*(R/r)^2*(5*(x(3)/r)^2-1));
dj_dt(5)=x(2)/x(1)*dj_dt(4);
dj_dt(6)=-(um*x(3)/r3)*(-(3/2)*J2*(R/r)^2*(5*(x(3)/r)^2-3));

return