function dr = n_bodies_function_Oc(x, t, mu, n)
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
a=zeros(n*3,1);
%Loop for the n bodies
for i=0:(n-1)
    for j=0:(n-1)
        if j~=i
	    dis = (((x(3*i+1)-x(3*j+1)))^2+((x(3*i+2)-x(3*j+2)))^2+((x(3*i+3)-x(3*j+3)))^2); 	     
	    for kk=1:3
              coo = (x(3*i+kk)-x(3*j+kk)); 	     
              a(3*i+kk)=a(3*i+kk) - mu(j+1) * coo / dis^(3/2) ;
            end
        end
    end
end
%        
dr=[x(3*n+1:end);a];
%
end
