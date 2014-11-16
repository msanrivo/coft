function tbd_pert = third_body_perturbations (et, pos, othr_params, Env)

% This function computes the gravity perturbations of selected third bodies
% from the knowledge of the S/C state vector in the MEE2000 reference
% frame and the ephemeris time.
% INPUTS:
%  - et     : Current ephemeris time [sec]
%  - pos    : Current S/C position vector in MEE2000 reference frame [km]
%  - othr_params:
%       .mu_tdb : Vector containing the third bodies gravitational
%                 potentials [km^3/s^2]
%  - Env    : Environment definition structure
% OUTPUTS:
%  - tbd_pert: Total third body acceleration in MEE2000 components [km/s^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------- START OF CODE -----------------------------
% Initialize third body perturbation vector
tbd_pert = [0;0;0];
% Central body code and string
cbody_str = Env.cbody;

for i = 1: Env.num_tdb
    
    % Gravitational potential of perturbing third body
    mu_tdb = othr_params.tdb_mu (i);
    
    % Retrieve the ephemerides of the perturbing celestial body in the
    % MEE2000 reference frame with respect to the central body
    tdb_str = Env.tdb_names{i}; % Third body string name
    [tdb_state,~] = cspice_spkezr(tdb_str, et, 'J2000', 'NONE', cbody_str);

    % Compute the thid body's gravity perturbation
    sc2tdb_vect    = tdb_state(1:3) - pos;
    sc2tdb_dist    = norm (sc2tdb_vect);
    cbody2tdb_vect = tdb_state(1:3);
    cbody2tdb_dist = norm (cbody2tdb_vect);
    pert           = mu_tdb * (sc2tdb_vect /(sc2tdb_dist^3) - ...
                               cbody2tdb_vect / (cbody2tdb_dist^3) );
                           
    % Sum up te perturbation of the ith perturbing body
    tbd_pert       = tbd_pert + pert;
end



