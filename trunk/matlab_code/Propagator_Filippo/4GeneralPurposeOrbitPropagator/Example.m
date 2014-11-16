% This script is an example of how to initialize SPICE data and prepare
% data for the propagator
close all;
clc;
clear all;

% Add the path to the SPICE (MICE) library and matlab source files
% THIS SHOULD BE THE PATH TO THE MICE FOLDER IN YOUR MACHINE!!!!
path2mice = 'C:\Users\fcichock\Dropbox\Astrodynamics\0reference_material\0code\5MicePackage';
% Load the required MICE kernels data
ME = load_spice_kernels ( path2mice );

% ALWAYS START YOUR SCRIPT WITH THE ABOVE LINES (PATH TO YOUR MICE FOLDER
% + CALL TO THE KERNELS LOAD FUNCTION)

% THE FOLLOWING IS AN EXAMPLE OF HOW TO PREPARE THE INPUT TO THE PROPAGATOR
%

% Initial epoch in ephemeris time
et0 = cspice_str2et( '1996 Dec 08 00:00:00' );
% Initial S/C state vector (MME2000)
mu_earth = cspice_bodvrd( 'EARTH', 'GM', 1 );
% Initial orbital elements (remember that these, for spice, are:
% E0(1) = pericentre radius
% E0(2) = eccentricity
% E0(3) = inclination
% E0(4) = right ascension of the ascending node
% E0(5) = argument of pericentre
% E0(6) = mean anomaly
% E0(7) = ephemeris time to which the above orbital elements refer to
% E0(8) = gravitational potential of the central body
E0 = [6878; 0.9 ; pi/4 ; 3/4*pi ; pi/2 ; pi ; et0 ; mu_earth];
S0      = zeros(7,1);
% Transformation from orbital elements to state vector using SPICE function
S0(1:6) = cspice_conics(E0, et0);
% Initial S/C mass (kg)
S0(7)   = 1000;

% Initialize some S/C parameters   
SC_params.cs_drag = 20.0;   % Cross section to drag [m2]
SC_params.cd      = 2.0;    % Drag coefficient
SC_params.cs_srp  = 20.0;   % Cross section to SRP [m2]
SC_params.rc_srp  = 1.25;   % SRP reflection coefficient
% Define the environment for the simulation
Env.cbody         = 'EARTH';  % Central body is the Earth (in this case)
% Perturbations to be considered
Env.insphe        = 1;      % Non spherical gravity flag (1:ON, 0:OFF)
Env.nsphe_deg     = 5;      % Non spherical gravity degree
Env.nsphe_ord     = 5;      % Non spherical gravity order
Env.idrag         = 2;      % Drag model flag (0:OFF, 1:simple model, 2:complex model)
Env.isrp          = 1;      % Solar radiation pressure flag (1:ON, 0:OFF)
Env.itdb          = 1;      % Third body perturbations flag (1:ON, 0:OFF)
Env.num_tdb       = 4;      % Number of third bodies for gravity perturbation
Env.tdb_names{1}  = 'SUN';
Env.tdb_names{2}  = 'MERCURY';
Env.tdb_names{3}  = 'VENUS';
Env.tdb_names{4}  = 'MARS';

% Propagation duration setting
deltat = 7*86400;
% Print out step for generating ephemeris
prnt_out_dt = 60*5;
% No stopping conditions other than elapsed time
stop_fun   = 'none';
% % Stop propagation at node crossing
% stop_fun = @(et,state)node_crossing_detection(et,state); 
% % Stop propagation at pericentre pass
% stop_fun = @(et, state)pericentre_pass_finder (et, state, mu_earth); 
% % Stop propagation at a distance of 10000 km from Earth's centre
% stop_fun = @(et, state)orbit_radius_crossing (et, state, 10000);

% Propagate with the Cowell's propagator
[SF, etf, states, times]  = cowell_propagator (S0, et0, deltat, ...
                                    SC_params, Env, prnt_out_dt, ...
									stop_fun);
  
% NEXT COMMENTED LINES SHOW HOW TO APPLY A MANOEUVRE AT THE FINAL EPOCH
% OF THE PROPAGATION
% % Apply a manoeuvre anti-parallel to the S/C inertial velocity at the end
% % of the previous propagation and of size 300 m/s
% DeltaV = SF(4:6)/norm(SF(4:6)) * 0.3;
% S0  = SF;
% et0 = etf;
% S0(4:6) = S0(4:6) - DeltaV;
% Isp = 300;
% g0  = 9.81;
% Update the S/C mass (due to the burn) with Tsiolkowsky equation
% S0(7) = S0(7)*exp(-DeltaV/(Isp*g0));
% % Propagate past the current stop condition (if necessary)
% [SF, etf, states_, times_] = cowell_propagator (S0, et0, 60, ...
%                                     SC_params, Env, 60, ...
% 									'none');
% S0 = SF;
% et0 = etf;
% [SF, etf, states2, times2]  = cowell_propagator (S0, et0, deltat, ...
%                                     SC_params, Env, prnt_out_dt, ...
% 									  stop_fun);

                                
                                
%% 
% BELOW IS AN EXAMPLE OF HOW TO PLOT THE RESULTS OF THE ABOVE SIMULATION
% Please observe that the output "states" is indeed a matrix with 7 rows
% and a number of columns equal to the number of printed time steps.
% If you are careful about the dimensions, Matlab allows you to "vectorize"
% function calls, as done in the following lines

% Reformulate dimensions of output matrices   
masses = states(7,:);
posvel = states(1:6,:);
% Convert the S/C state vectors into a S/C elements vectors (a set of 
% "nsteps" state vectors of 6 components is transformed into a set of 
% "nsteps" orbital elements vectors of 8 components)
sc_elems = cspice_oscelt( posvel, times, mu_earth );
                                                    
% 3-D PLOT OF THE ORBIT
figure(1);
dim = size(states,2);
[X,Y,Z] = ellipsoid(0,0,0,6378,6378,6378);
h = surfl(X, Y, Z);
hold on;
set (h, 'FaceColor', [0.7 0.7 0.7]);
hold on;
plot3(states(1,:),states(2,:),states(3,:),'r');
% hold on;
% plot3(states2(1,:),states2(2,:),states2(3,:),'--k');
hold on;
plot3([0,6378*2], [0,0], [0,0], 'k','LineWidth',2);
hold on;
plot3([0,0], [0,6378*2], [0,0], 'k','LineWidth',2);
hold on;
plot3([0,0], [0,0], [0,6378*2], 'k','LineWidth',2);
hold on;
axis equal;
xlabel ('X coordinate of the MEE2000 reference frame [km]','Fontsize',24);
ylabel ('Y coordinate of the MEE2000 reference frame [km]','Fontsize',24);
zlabel ('Z coordinate of the MEE2000 reference frame [km]','Fontsize',24);
set(gca,'Fontsize',24);

% Compute the analytical change of RAAN due to the J2 effect
J2 = 1.08262668355*1e-3;
RAAN0 = sc_elems(4,1);
Delta_RAAN = sc_elems(4,end)-RAAN0;
p = sc_elems(1,1)*(1+sc_elems(2,1));
sma = sc_elems(1,1)/(1-sc_elems(2,1));
inc = sc_elems(3,1);
N_orbits   = (times - et0) ./ (2*pi*sqrt(sma^3/mu_earth));
RAAN_analytical = RAAN0 - N_orbits * 3 * pi * J2 * (6378/p)^2 * cos (inc);
Delta_RAAN_analytical = RAAN_analytical (end) - RAAN0;

% EVOLUTION OF THE RIGHT ASCENSION OF THE ASCENDING NODE
figure(2);
plot((times-et0)/3600, sc_elems(4,:)*180/pi, 'r', 'LineWidth', 2);
hold on;
plot((times-et0)/3600, RAAN_analytical*180/pi, '--k', 'LineWidth', 2); 
xlabel ('Hours since the start of the simulation','Fontsize',24);
ylabel ('Right Ascension of the ascending node [deg]','Fontsize',24);
set(gca,'Fontsize',24);

% EVOLUTION OF THE OSCULATING SEMI-MAJOR AXIS
figure(3);
plot((times-et0)/3600, sc_elems(1,:) ./ (1-sc_elems(2,:)), 'r');
xlabel ('Hours since the start of the simulation','Fontsize',24);
ylabel ('Semi-major axis [km]','Fontsize',24);
set(gca,'Fontsize',24);

% EVOLUTION OF THE OSCULATING ECCENTRICITY 
figure(4);
plot((times-et0)/3600, sc_elems(2,:), 'r');
xlabel ('Hours since the start of the simulation','Fontsize',24);
ylabel ('Eccentricity','Fontsize',24);
set(gca,'Fontsize',24);

% EVOLUTION OF THE OSCULATING INCLINATION
figure(5);
plot((times-et0)/3600, sc_elems(3,:)*180/pi, 'r');
xlabel ('Hours since the start of the simulation','Fontsize',24);
ylabel ('Inclination with respect to Earth equator [deg]','Fontsize',24);
set(gca,'Fontsize',24);

% EVOLUTION OF THE OSCULATING ARGUMENT OF PERICENTRE
figure(5);
plot((times-et0)/3600, sc_elems(5,:)*180/pi, 'r');
xlabel ('Hours since the start of the simulation','Fontsize',24);
ylabel ('Argument of pericentre [deg]','Fontsize',24);
set(gca,'Fontsize',24);

% Show the difference between the theoretical RAAN drift due to J2 and 
% the real one
( Delta_RAAN - Delta_RAAN_analytical ) / Delta_RAAN
