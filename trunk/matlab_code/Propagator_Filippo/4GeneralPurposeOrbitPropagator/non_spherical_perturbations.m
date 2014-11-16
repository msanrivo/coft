function non_sphere_pert = non_spherical_perturbations ...
                           (et, pos, othr_params, Env)

% This function computes the non-spherical gravity perturbations of the 
% Earth in the MEE2000 reference frame, given the current S/C state vector
% and the ephemeris time
% INPUTS:
%  - et     : Current ephemeris time [sec]
%  - pos    : Current S/C position vector in MEE2000 reference frame [km]
%  - othr_params:
%       .Re_cen: Equatorial radius of the central body [km]
%       .mu_cen: Gravitational potential of the centralbody [km^3/s^2]
% OUTPUTS:
%  - non_sphere_pert: Non-spherical gravity acceleration in MEE2000 components 
%                     [km/s^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earth radius parameter [km]
Mu  = othr_params.mu_cen;      % [km^3/s^2]
Re  = othr_params.Re_cen;      % [km]

% Rotate the S/C state vector to the rotating equatorial reference frame of
% the central body 
cbody_rot_frame = ['IAU_',Env.cbody];
try
    rot_MEE2000_2_ROTFRAME = cspice_pxform('J2000', cbody_rot_frame, et);
catch
    rot_MEE2000_2_ROTFRAME = [1,0,0;0,1,0;0,0,1];
end
    
% Compute S/C position vector in such rotating reference frame
pos_fixed = rot_MEE2000_2_ROTFRAME * pos;
x = pos_fixed(1);
y = pos_fixed(2);
z = pos_fixed(3);
r = norm(pos_fixed);

lat       = asin(z/r);
long      = atan2(y, x);
% Acceleration function takes inputs in [m], so we should pass the inputs
% accordingly!! We want output in [km/s] though!!!
degree     = Env.nsphe_deg;
order      = Env.nsphe_ord;
pert_fixed = non_sphere_model (othr_params.C, othr_params.S, lat, long, ...
                              r, Re, Mu, degree, order);

% Rotate perturbative acceleration from body fixed components to MEE2000 
% components (observe that it is equivalent to pre-multiply for the 
% transverse matrix of the MEE2000 to rotating frame transformation
non_sphere_pert = rot_MEE2000_2_ROTFRAME' * pert_fixed;

