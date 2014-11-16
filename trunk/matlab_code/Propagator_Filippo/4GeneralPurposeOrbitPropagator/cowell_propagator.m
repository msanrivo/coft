% ------------------- NUMERICAL INTEGRATION FUNCTION ----------------------
% This function propagates the state of a S/C in the non-rotating Mean 
% Earth equator reference frame (SPICE "J2000" reference frame). The
% required SPICE data kernels must be loaded before calling this function
% using the "load_spice_kernels.m" function provided.
% You can propagate the orbit in any inertial reference frame, provided
% that the perturbations considered do not depend on the attitude of the
% central body. In that case, the code assumes by default that the S/C
% state vector is given in MEE2000 (SPICE "J2000") components
%
% -------------------------------------------------------------------------
% INPUTS:
%  - S0 [1:7] : Initial spacecraft state vector in [km] and [km/s] in the
%               MEE2000 reference frame (Mean Earth Equator at 2000 Jan 1st 
%               12:00:00). Last component is the initial S/C mass [kg]
%
%  - SC_params: Structure containing various info on the S/C properties
%       .cs_drag : Drag cross section [m^2]
%       .cd      : Drag coefficienct
%       .cs_srp  : SRP cross section [m^2]
%       .rc_srp  : SRP reflective coefficient [0, 2]
%
%  - Env      : Structure containing information on the environment to
%               be simulated:
%      .cbody    : Upper case name of the central body of the propagation:
%                  - 'SUN', 'MERCURY', 'VENUS', 'EARTH', 'MOON', 'MARS',
%                    'JUPITER','SATURN','NEPTUNE','URANUS' etc...
%                  - 'USER_DEFINED': User defined central body
%      .cbody_mu : If Env.cbody = 'USER_DEFINED', then this represents the 
%                  gravitational potential of the central body of interest 
%                  [km^3/s^2]
%      .isrp     : Flag to activate solar radiation pressure perturbation
%                  - 0 : No SRP considered
%                  - 1 : SRP considered with Cannonball model
%      .insphe   : Flag to activate non spherical gravity model
%                  - 0 : No effects of non-spherical gravity
%      .insphe_deg: In case of .insphe = 1, this is the maximum degree of
%                  the gravity coefficients to be considered
%      .insphe_ord: In case of .insphe = 1, this is the maximum order of 
%                  gravity coefficients to be considered
%      .idrag    : Flag to activate the drag perturbation effects
%                  - 0 : No effects of drag
%                  - 1 : Simplified drag model
%                  - 2 : More complex drag model
%      .itdb     : Flag to activate the third body gravity perturbations
%      .num_tdb  : Number of third bodies to be included (max is 10)
%      .tbd_names: Vector containing the SPICE string names of the 
%                  celestial bodies to be included as third bodies
%
%  - et0      : Initial ephemeris time [sec]. Elapsed seconds since TT
%               2000 Jan 1st 12:00:00. This is the american convention for
%               time counting in mission analysis problems
%
%  - deltat   : Propagation duration [sec]
%
%  - prnt_out_dt: Print-out step for the S/C state vector
%
%  - stop_fun : Function handle to the function determining the exit
%               conditions (if different from the total propagation time).
%               This input is optional !!!
%
%
% OUTPUTS:
%  - SF [1:7] : Final spacecraft state vector in [km] and [km/s] in the
%               MEE2000 reference frame (ECI at 2000 Jan 1st 12:00:00).
%               Last component is the final S/C mass [kg]
%  - etf      : Final ephemeris time corresponding to SF [s]
%  - states [1:7,dim]: S/C state vectors at the print out steps. The 
%               number of columns is equal to the number of time steps
%               that have been processed
%  - times [dim]: Ephemeris times row vector [s]. Its size is equal to
%               the number of time steps that have been processed
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SF, etf, states, times] = cowell_propagator (S0, et0, deltat, ...
                                    SC_params, Env, prnt_out_dt, stop_fun)

% Initial S/C state vector in MEE2000 ("J2000" reference frame)
state = S0;
% Final ephemeris time
etf = et0 + deltat;

% Initialize parameters structure for the use of the derivatives function
% Required planetary constants
if ( strcmp(Env.cbody,'USER_DEFINED')== 0)
    othr_params.mu_cen     = cspice_bodvrd( Env.cbody, 'GM', 1 );
    radii                  = cspice_bodvrd( Env.cbody, 'RADII', 3 );
    othr_params.Re_cen     = radii(1);
else
    othr_params.mu_cen     = Env.cbody_mu;
end

% Initialize the sidereal rotation period of the central body for the
% atmospheric drag perturbation and send warning messages, if required
warning1 = ['Atmosphere model or sidereal rotation period is not ', ...
            'available for this central body'];
if (strcmp(Env.cbody,'EARTH'))
    othr_params.cbody_spin = 7.292115*10^-5; % Earth spin rate [Rad/s]
elseif ( strcmp(Env.cbody,'MERCURY')) 
    disp (warning1);
    sid_rot_period = 58.6462*86400;
    othr_params.cbody_spin = 2*pi/sid_rot_period; % Mercury spin rate [Rad/s]
elseif ( strcmp(Env.cbody,'VENUS'))  
    disp (warning1);   
    sid_rot_period = 243.0187*86400;
    othr_params.cbody_spin = 2*pi/sid_rot_period; % Venus spin rate [Rad/s]
elseif ( strcmp(Env.cbody,'MARS'))
    disp (warning1); 
    sid_rot_period = 1.02595675*86400;
    othr_params.cbody_spin = 2*pi/sid_rot_period; % Mars spin rate [Rad/s]
elseif ( strcmp(Env.cbody,'JUPITER'))
    disp (warning1);  
    sid_rot_period = 0.41007*86400;
    othr_params.cbody_spin = 2*pi/sid_rot_period; % Jupiter spin rate [Rad/s]
elseif ( strcmp(Env.cbody,'MOON'))
    disp (warning1);
    sid_rot_period = 27.321661*86400;
    othr_params.cbody_spin = 2*pi/sid_rot_period; % Moon spin rate [Rad/s]
else
    disp (warning1);
end
    
% Initialize the gravitational parameters of required third bodies
for i = 1:Env.num_tdb
    othr_params.tdb_mu(i) = cspice_bodvrd( Env.tdb_names{i}, 'GM', 1);
end

% Physical constants of interest
% Solar radiation pressure [Pa] at a distance from the Sun of 1 AU   
othr_params.pre_srp0    = 4.57*10^-6; 

% Gravitational coefficients of the central body gravity
if ( strcmp(Env.cbody,'EARTH'))
	[othr_params.C, othr_params.S] = ini_coeff_earth;
else
    disp('No gravity coefficients are available for this central body');
end

% Function that computes the inertial acceleration components on the S/C
tot_derivs = @(et,state)state_derivs(et, state, SC_params, ...
                                     othr_params, Env);

% Options of the Runge Kutta solver
% Maximum integration time step is 1/50 of the orbital period
sc_elems = cspice_oscelt( state(1:6), et0, othr_params.mu_cen );
sma      = sc_elems(1) ./ (1-sc_elems(2));
orb_per  = 2*pi*sqrt(sma^3/othr_params.mu_cen);
max_step = 1/50 * orb_per;

% Additional exit condition of the Cowell's propagator (if it is reached
% before the final propagation time)
options = odeset('RelTol',1e-9,'AbsTol',1e-9, 'MaxStep',max_step);
% Set events function to the input function handle
if ( strcmp(stop_fun, 'none') == 0 )
    options.Events = stop_fun;
end
% Avoid errors when the print out step has been set higher than the
% propagation duration
if (prnt_out_dt > etf-et0)
    prnt_out_dt = etf-et0;
end
% ---------- SOLVE FOR THE TRAJECTORY WITH AN ODE45 INTEGRATOR ------------
[times,states] = ode45(tot_derivs,et0:prnt_out_dt:etf,state,options);       

% Reformulate output vector dimensions
times  = times';
states = states';

% Update the final S/C state value and ephemeris time
SF = states(1:7,end);
etf = times(end);

return;
