% Calcula el efecto perturbador de la presion solar
%
% function df_dt = flux( x ) 
% 
% Author:   
% Date:     
 
% Entrada:
%       
%	x, vector posición/velocidad
      
% Salida:
%
%       df_dt, vector con la v_pert y la a_pert
%
%           df_dt(1), v_pert(1)
%           df_dt(2), v_pert(2)
%           df_dt(3), v_pert(3)
%           df_dt(4), a_pert(1)
%           df_dt(5), a_pert(2)
%           df_dt(6), a_pert(3)

function [ df_dt ] = flux( x )

jd = jday(2005,07,21,  09,49,08);

[rsun,time,decl] = sun ( jd );

x_sat(1)=rsun(1)-x(1);

x_sat(2)=rsun(2)-x(2);

x_sat(3)=rsun(3)-x(3);

r_sat=norm(x_sat(1:3));

a_sr=3.41e-7; %-4.51e-7; %[m/s^2] a_sr=p_sr*c_r*A/m;

df_dt(1)=0;
df_dt(2)=0;
df_dt(3)=0;
df_dt(4)=-a_sr*x_sat(1)/r_sat;
df_dt(5)=-a_sr*x_sat(2)/r_sat;
df_dt(6)=-a_sr*x_sat(3)/r_sat;

return