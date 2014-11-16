function [derivs] = state_derivs (et, state, SC_params, othr_params, Env)

% This function computes the inertial acceleration acting at time "et" on 
% the S/C in [km/s^2] and the rate of change of the S/C mass as a function
% of the current ephemeris time, S/C state and input parameters
% INPUTS:
%  - et     : Current ephemeris time [sec]
%  - state  : Current S/C state in MEE2000 reference frame [km, km/s]
%  - SC_params: Structure containing a number of S/C parameters
%  - othr_params: Structure containing a number of environmental parameters
%  - Env    : Environment definition parameters
% OUTPUTS:
%  - derivs: Total derivatives of the S/C state vector
%     - derivs (1:3) : Inertial position derivatives (MEE2000)
%     - derivs (4:6) : Inertial velocity derivatives (MEE2000)
%     - derivs (7)   : S/C mass derivative [kg/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize derivative vector (7 components to include mass rate as well)
derivs = zeros(7,1);

% Compute central gravity component
mu_cen = othr_params.mu_cen;
radius = norm (state(1:3));
central_gravity = - mu_cen * state(1:3) / (radius^3);

% Compute gravity perturbations due to zonal terms up to degree 6
if (Env.insphe == 1)
    non_sphere_pert = non_spherical_perturbations (et,state(1:3),...
                                                   othr_params, Env);
else
    non_sphere_pert = [0;0;0];
end

% Compute third body gravity perturbations
if (Env.itdb == 1) 
    tdb_pert = third_body_perturbations (et, state(1:3), othr_params, Env);
else
    tdb_pert = [0;0;0];
end 

% Compute solar radiation pressure effects
if (Env.isrp == 1)
    srp_pert = srp_perturbations (et, state, SC_params, othr_params, Env);
else
    srp_pert = [0;0;0];
end

% Compute drag perturbation acceleration
if (Env.idrag > 0)
    drag_pert = drag_perturbations(et, state, SC_params, othr_params, Env);
else
    drag_pert = [0;0;0];
end

% Total position derivatives
derivs(1:3) = state(4:6);
% Total velocity derivatives
derivs(4:6) = central_gravity + non_sphere_pert + tdb_pert + srp_pert + ...
              drag_pert;
% S/C mass derivatives     
derivs(7)   = 0.0;
