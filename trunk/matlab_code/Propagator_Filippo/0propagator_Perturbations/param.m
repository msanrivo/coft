% Calcula la variacion que se produce en los parametros orbitales.
%
% function dx_dt = param( t, x, um, p0, R, perturb_J2, perturb_drag, perturb_flux ) 
% 
% Autor: Paloma Belen Mangas Carbajo   
  
 
% Entrada:
%       
%	x, parametros orbitales
%   um, parametro gravitacional
%   p0, parametros orbitales iniciales
%   R,radio de la Tierra
%
% Salida:
%
%       dx_dt, vector con las variaciones de los parametros
%
%	dx_dt(1)=da_dt
%	dx_dt(2)=de_dt  
%	dx_dt(3)=di_dt 
%	dx_dt(4)=dO_dt	
%	dx_dt(5)=dw_dt 
%	dx_dt(6)=dM_dt 

function [ dx_dt ]=param( t, x, um, p0, R, perturb_J2, perturb_drag )

n=sqrt(um/x(1)^3);

[ dj2p_dt ]=J2_param( x, um, p0, R, n );
[ ddp_dt ]=drag_param( x, um );

if perturb_J2 == 1
[ dj2p_dt ] = J2_param( x, um, R, r, r3 );
else
    dj2p_dt = zeros([6,1]); 
end
if perturb_drag == 1
[ ddp_dt ] = drag_param( x, r );
else
    ddp_dt = zeros([6,1]); 
end

dx_dt(1)=perturb_J2*dj2p_dt(1)+perturb_drag*ddp_dt(1);
dx_dt(2)=perturb_J2*dj2p_dt(2)+perturb_drag*ddp_dt(2);
dx_dt(3)=perturb_J2*dj2p_dt(3)+perturb_drag*ddp_dt(3);
dx_dt(4)=perturb_J2*dj2p_dt(4)+perturb_drag*ddp_dt(4); 
dx_dt(5)=perturb_J2*dj2p_dt(5)+perturb_drag*ddp_dt(5); 
dx_dt(6)=(180/pi)*n+perturb_J2*dj2p_dt(6)+perturb_drag*ddp_dt(6); 

dx_dt=dx_dt';

return