function dr = n_bodies_function(t,x,mu,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function integrates the position and velocity
%
% Author: Pablo A. Machuca Varela
%
% Input:
% t              [1x1]     Temporal value
% x              [6*nx1]   [x; y; z; vx; vy; vz] for the value t in time
% mu             [nx1]     Gravitational parameters
% n              [1x1]     Number of bodies
%
% Output:
% a         [3*nx1]      [ax; ay; az]                Acceleration
% dr        [6*nx1]      [dx; dy; dz; dvx; dvy; dvz] Derivatives of position and velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Preallocation
a = zeros(n*3,1);
%
%Loop for the n bodies
for i=0:(n-1)
    for j=0:(n-1)
        if j~=i
        a(3*i+1) = a(3*i+1) - (mu(j+1)*(x(3*i+1)-x(3*j+1)))/(((x(3*i+1)-x(3*j+1)))^2+((x(3*i+2)-x(3*j+2)))^2+((x(3*i+3)-x(3*j+3)))^2)^(3/2); %ax = mu*(x2-x1)/(r2-r1)^(3/2)
        a(3*i+2) = a(3*i+2) - (mu(j+1)*(x(3*i+2)-x(3*j+2)))/(((x(3*i+1)-x(3*j+1)))^2+((x(3*i+2)-x(3*j+2)))^2+((x(3*i+3)-x(3*j+3)))^2)^(3/2); %ay = mu*(y2-y1)/(r2-r1)^(3/2)
        a(3*i+3) = a(3*i+3) - (mu(j+1)*(x(3*i+3)-x(3*j+3)))/(((x(3*i+1)-x(3*j+1)))^2+((x(3*i+2)-x(3*j+2)))^2+((x(3*i+3)-x(3*j+3)))^2)^(3/2); %az = mu*(z2-z1)/(r2-r1)^(3/2)
        end 
    end
end
%      
dr = [x(3*n+1:end); a];
%
end
