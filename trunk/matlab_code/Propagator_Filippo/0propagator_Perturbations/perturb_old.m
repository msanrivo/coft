%Calcula la variacion que se produce en r y v
%
% function dx_dt = perturb( t, x, um, r0, perturb_J2, perturb_drag ) 
% 
% Autor: Paloma Belen Mangas Carbajo      
 
% Entrada:
%       
%	x, vector posición/velocidad
%   um, parametro gravitacional
%	r0, sqrt(x0(1)^2+x0(2)^2+x0(3)^2)
%       		   
% Salida:
%
%   dx_dt, vector con v y a
%
%	dx_dt(1)=v(1)
%	dx_dt(2)=v(2)
%	dx_dt(3)=v(3) 
%	dx_dt(4)=a(1)
%	dx_dt(5)=a(2) 
%	dx_dt(6)=a(3) 


function [ dx_dt ] = perturb( t, x, um, perturb_J2, perturb_drag, perturb_flux )

R=6378145;

r3=(norm(x(1:3)))^3;
r=(norm(x(1:3)));

if perturb_J2 == 1
[ dj_dt ] = J2( x, um, R, r, r3 );
else
    dj_dt = zeros([6,1]); 
end
if perturb_drag == 1
[ dd_dt ] = drag( x, r );
else
    dd_dt = zeros([6,1]); 
end
if perturb_flux == 1
[ df_dt ] = flux( x );
else
    df_dt = zeros([6,1]); 
end

dx_dt(1)=x(4)+dj_dt(1)*perturb_J2+dd_dt(1)*perturb_drag+perturb_flux*df_dt(1);
dx_dt(2)=x(5)+dj_dt(2)*perturb_J2+dd_dt(2)*perturb_drag+perturb_flux*df_dt(2);
dx_dt(3)=x(6)+dj_dt(3)*perturb_J2+dd_dt(3)*perturb_drag+perturb_flux*df_dt(3);
dx_dt(4)=-um*x(1)/r3+dj_dt(4)*perturb_J2+dd_dt(4)*perturb_drag+perturb_flux*df_dt(4);
dx_dt(5)=-um*x(2)/r3+dj_dt(5)*perturb_J2+dd_dt(5)*perturb_drag+perturb_flux*df_dt(5);
dx_dt(6)=-um*x(3)/r3+dj_dt(6)*perturb_J2+dd_dt(6)*perturb_drag+perturb_flux*df_dt(6);

dx_dt=dx_dt';

return
