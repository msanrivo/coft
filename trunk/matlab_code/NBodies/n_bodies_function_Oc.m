function dr = n_bodies_function_Oc(x, t, mu, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION n_bodies_function
% 
% AUTHOR:    Pablo Machuca
% DATE:      15/02/2013
% VERSION:   1.0
%
% This function computes the derivatives of the state 
% of n bodies ...
%
% Input:
%       
%	x  [nx6] vector of states...
%       mu [nx1] array of gravitational parameters...
%       ....
%	   
% Output:
%
%       dr [nx6] vector of derivatives
%
% Usage: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=zeros(n*3,1);
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
