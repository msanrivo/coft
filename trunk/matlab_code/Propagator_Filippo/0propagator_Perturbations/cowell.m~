% Integrador de r, v a partir de r0 y v0
%
% function cowell 
% 
% Autor: Paloma Belen Mangas Carbajo   
%
% 
% Entrada:
%       
%   x0, vector posición/velocidad
%   q, tiempo de simulacion / um 
%   perturb_J2, activa/desactiva esta perturbacion
%   perturb_drag, activa/desactiva esta perturbacion
%   perturb_flux, activa/desactiva esta perturbacion
%	   
% Salida:
%
%   x, vector con la x_pert y la v_pert
%
%           x(1), x_pert(1)
%           x(2), x_pert(2)
%           x(3), x_pert(3)
%           x(4), v_pert(1)
%           x(5), v_pert(2)
%           x(6), v_pert(3)

format long e 

[x0,q,perturb_J2,perturb_drag,perturb_flux] = um_tf_xo_cowell;

um=q(1)*q(2);
tf=q(3);

tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
tspan=(0:10:tf);
[t,x]=ode45( @perturb,tspan,x0,options,um,perturb_J2,perturb_drag,perturb_flux );

z0=x;

for i=1:length(t)

q0=z0(i,:);

p0(i,:)=rv2coe(q0);

z1(i,:)=z0(i,1:3);
z2(i,:)=z0(i,4:6);
r(i)=norm(z0(i,1:3));
v(i)=norm(z0(i,4:6));
u(i)=um/r(i);
e(i)=v(i).^2/2 - u(i);
h(i,:)=cross(z1(i,:),z2(i,:));
hk(i)=norm(h(i,1:3));

end

fi = fopen('salida_cowell.txt', 'w');

fprintf(fi, '%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n', p0'); 

fclose(fi);

r_x=x(end,1)
r_y=x(end,2)
r_z=x(end,3)
v_x=x(end,4)
v_y=x(end,5)
v_z=x(end,6)
 

%Plots

figure(1)

%hk=componente perpendicular del momento angular h=rxv
%(hk - hk_0)/hk_0%

semilogy(t, abs(hk-hk(1))/abs(hk(1)) , 'ko-.')
grid on
hold on
xlabel('t (s)');
ylabel('(hk - hk_0)/hk_0 ');

figure(2)

%e=energía%
%(e - e_0)/e_0%

semilogy(t, abs(e-e(1))/abs(e(1)) , 'ko-.')
grid on
hold on
xlabel('t (s)');
ylabel('(E - E_0)/E_0');