% This function computes the non-spherical gravity perturbations up to
% degree N
function [pert_non_sphere] = non_sphere_model (C, S, lat, long, r, ...
                             Re, Mu, degree, order)

% Compute the derivative of the radius, latitude and longitude with respect
% to the radius vector
x = abs(r*cos(lat))*cos(long);
y = abs(r*cos(lat))*sin(long);
z = r*sin(lat);

% Derivatives with respect to a vector are row vectors
% Radial distance derivatives with respect to the cartesian coordinates
dr_drvec = [1/r * x , 1/r * y, 1/r * z];
% Latitude derivatives with respect to the cartesian coordinates
dlat_drvec(1,1) = 1/sqrt(x^2+y^2)* (-x*z/r^2);
dlat_drvec(1,2) = 1/sqrt(x^2+y^2)* (-y*z/r^2);
dlat_drvec(1,3) = 1/sqrt(x^2+y^2)* (-z*z/r^2 + 1);
% Longitude derivatives with respect to the cartesian coordinates
dlong_drvec(1,1) = 1/(x^2+y^2)* ( - y );
dlong_drvec(1,2) = 1/(x^2+y^2)* ( + x );
dlong_drvec(1,3) = 0;


% Derivative of the non-spherical part of the potential function with 
% respect to the radial distance

P = zeros(6,6);
gamma  = sin (lat);
P(1,1) = 1; 
P(2,1) = gamma;
P(2,2) = cos(lat);

dU_dr = 0;
dU_dlat = 0;
dU_dlong = 0;
for l = 3 : degree+1
    for m = 1: min(order+1, l)
        % Compute Legendre Function
        P(l,m)   = legendre(P,l,m,sin(lat));
        P(l,m+1) = legendre(P,l,m+1,sin(lat));
        
        % Coefficients as they appear in the formulas
        ll = l - 1;
        mm = m - 1;
        % Increase in the potential function derivative with respect to the
        % radius
        ddU_dr = - Mu/(r^2) * (Re/r)^(ll) * (ll+1) * P(l,m) * ...
                  (C(l,m)*cos(mm*long) + S(l,m) * sin(mm*long));
              
        % Increase in the potential function derivative with respect to the
        % latitude
        ddU_dlat = + Mu/r * (Re/r)^(ll) * ...
                   ( P(l,m+1) - mm * tan(lat) * P(l,m) ) * ...
                   (C(l,m) * cos( mm*long) + S(l,m)*sin( mm*long));
        
        % Increase in the potential function derivative with respect to the
        % longitude
        ddU_dlong = + Mu/r * (Re/r)^(ll)* mm * P(l,m) * ...
                   (S(l,m) * cos( mm*long) + C(l,m)*sin( mm*long));       
            
        % Update derivatives of the potential with newly computed terms       
        dU_dr    = dU_dr + ddU_dr;
        dU_dlat  = dU_dlat + ddU_dlat;
        dU_dlong = dU_dlong + ddU_dlong;
        
               
    end
end

pert_non_sphere = dU_dr * dr_drvec' + ...
                  dU_dlat * dlat_drvec' + ...
                  dU_dlong * dlong_drvec';





