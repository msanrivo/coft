% Pericentre pass finder (only for elliptical orbits)
function [lookfor stop direction] = pericentre_pass_finder ...
                                    (et, state, mu_central_body)

% Orbit radius
r = sqrt(state(1)^2+state(2)^2+state(3)^2);
% Compute the osculating orbital elements
elems = cspice_oscelt( state(1:6), et, mu_central_body );
% eccentriciy of the orbit
e = elems(2);
% Semi-major axis of the orbit
SMA = elems(1)/(1-e); 
% Semi-latum rectum
p = SMA*(1-e^2);
% Mean anomaly
M   = elems(6);
% Eccentric anomaly is then given by:
arg = (p/r-1)/e;
if ( abs(arg) >= 1.0)
    arg = 1*sign(arg);
end
if (M <= pi)
    TRA = acos ( arg );
else
    TRA = 2*pi - acos ( arg );
end
% Always consider true anomaly between -pi and pi (to have a continuous
% transition at the pericentre)
if (TRA > pi)
    TRA = TRA - 2*pi;
end

% Stopping criterion
lookfor   = TRA ;  %Searches for this expression set to 0 (pericentre pass)
stop      = 1;   %Stop when event is located
direction = 1;   %Specifiy direction of motion at event (from negative E)