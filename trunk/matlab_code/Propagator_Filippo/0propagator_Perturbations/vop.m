% Integrador de parametros_orbitales a partir de parametros_orbitales0
%
% function vop
% 
% Autor: Paloma Belen Mangas Carbajo  
 
% Entrada:
%       
%   p0, vector parametros_orbitales0
%   q, tiempo de simulacion / um 
%   perturb_J2, activa/desactiva esta perturbacion
%   perturb_drag, activa/desactiva esta perturbacion
%       	   
% Salida:
%
%       x, vector parametros_orbitales0 (a, e, i, omega, w, M)
%
%           x(1), a
%           x(2), e
%           x(3), i
%           x(4), omega
%           x(5), w
%           x(6), M

format long e 

[p0,q,perturb_J2,perturb_drag] = um_tf_xo;

um=q(1)*q(2);
tf=q(3);
x0=p0;
R=6378145;

tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
tspan=(0:10:tf);

[t,x]=ode45(@param,tspan,x0,options,um,p0,R,perturb_J2,perturb_drag);

a=x(end,1)
e=x(end,2)
i=x(end,3)
O=((x(end,4)/360)-floor(x(end,4)/360))*360
w=((x(end,5)/360)-floor(x(end,5)/360))*360
M=((x(end,6)/360)-floor(x(end,6)/360))*360

for i=1:length(t)

q0(i,1)=x(i,1);
q0(i,2)=x(i,2);
q0(i,3)=x(i,3);
q0(i,4)=((x(i,4)/360)-floor(x(i,4)/360))*360;
q0(i,5)=((x(i,5)/360)-floor(x(i,5)/360))*360;
q0(i,6)=((x(i,6)/360)-floor(x(i,6)/360))*360;

end

fi = fopen('salida.txt', 'w');

fprintf(fi, '%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n', q0'); 

fclose(fi);