% Calcula la variacion que se produce en los parametros orbitales debido a J2.
%
% function dj2p_dt = J2_param( x, um, p0, R, n )                                                                                                                                             
% 
% Autor: Paloma Belen Mangas Carbajo  
    
 
% Entrada:
%       
%       x, parametros orbitales
%       um, parametro gravitacional
%       p0, parametros orbitales iniciales
%	    R, radio medio de la Tierra
%       n, revoluciones
%
%	   
% Salida:
%
%       dj2p_dt, vector con las variaciones de los parametros
%
%	dx_dt(1)=da_dt
%	dx_dt(2)=de_dt  
%	dx_dt(3)=di_dt 
%	dx_dt(4)=dO_dt	
%	dx_dt(5)=dw_dt 
%	dx_dt(6)=dM_dt 

function [ dj2p_dt ]=J2_param( x, um, p0, R, n )

J2=0.00108248;
p=x(1)*(1-(x(2)^2));

dj2p_dt(1)=0; %(p0(1)*3*J2*(R^2)*(sqrt(1-x(2)^2)*(3*sind(x(3))^2-2)))/(4*p^2); 
dj2p_dt(2)=0; %((2*(1-x(2)))/(3*n))*n_dot;
dj2p_dt(3)=0;
dj2p_dt(4)=-((3*n*(180/pi)*R^2*J2*cosd(x(3)))/(2*p^2)); %Vallado pag 647
dj2p_dt(5)=(3*n*(180/pi)*R^2*J2*(4-5*sind(x(3))^2))/(4*p^2); %Vallado pag 647
dj2p_dt(6)=-(3*n*(180/pi)*J2*(R^2)*(sqrt(1-x(2)^2)*(3*sind(x(3))^2-2)))/(4*p^2); %Vallado pag 647

return