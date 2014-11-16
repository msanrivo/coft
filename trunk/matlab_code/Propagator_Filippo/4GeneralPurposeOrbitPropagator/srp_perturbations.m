function srp_pert = srp_perturbations (et, state, SC_params, ...
                                       othr_params, Env)

% This function computes the solar radiation perturbation, considering a
% simplified cannonball model.
% INPUTS:
%  - et     : Current ephemeris time [sec]
%  - pos    : Current S/C position vector in MEE2000 reference frame [km]
%  - SC_params:
%       .cs_srp : S/C cross section to solar radiation [m^2]
%       .rc_srp : S/C reflective coefficient [0:2]
%  - othr_params: Structure containing physical parameters
%  - Env    : Structure containing environment definition flags
% OUTPUTS:
%  - srp_pert: Vectorial acceleration due to SRP in MEE2000 components
%             [km/s^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solar radiation pressure at 1 AU
srp_pressure_1AU = othr_params.pre_srp0;
AU = 149597871;  % Astronomical unit in km
% Reflective coefficient for SRP
rc_srp = SC_params.rc_srp;
% Cross section to SRP
cs_srp = SC_params.cs_srp;

pos = state(1:3);
mass = state(7);

% Compute the distanceof the S/C from the Sun
cbody_str   = Env.cbody;
[sun_state,  ~] = cspice_spkezr('Sun', et, 'J2000', 'NONE', cbody_str);
sun2sc_vect = pos - sun_state(1:3);
sun2sc_dist = norm(sun2sc_vect);
sun2sc_vers = sun2sc_vect/sun2sc_dist;

% Compute effective SRP pressure at the real S/C distance from the Sun
srp_pressure = srp_pressure_1AU * ( (AU/sun2sc_dist)^2 );

% Compute vector perturbative acceleration due to SRP [in km/s^2]
srp_pert = - sun2sc_vers * srp_pressure * rc_srp * cs_srp / mass *1e-3;
