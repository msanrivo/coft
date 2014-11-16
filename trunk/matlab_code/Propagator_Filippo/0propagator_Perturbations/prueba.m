function [f, tf, x0, perturb_J2, perturb_drag, perturb_flux]=prueba

%f=1    formato TLE inputs por pantalla
%f=2    formato ECI
%f=3    FORMATO COE inputs por pantalla

f=2; 

tf=86400;

x0=[-2436450 -2436450 6891037.9 5088.611 -5088.611 0];

perturb_J2 = 0;
perturb_drag = 0;
perturb_flux = 0;