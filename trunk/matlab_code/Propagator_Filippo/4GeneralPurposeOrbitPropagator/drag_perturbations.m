function drag_pert = drag_perturbations (et, state, SC_params, ...
                                          othr_params, Env)

% This function computes the drag perturbation, considering a relatively
% complex atmospheric model and the true velocity wrt the atmosphere
% INPUTS:
%  - et     : Current ephemeris time [sec]
%  - state  : Current S/C state in MEE2000 reference frame [km, km/s, kg]
%  - SC_params:
%       .cs_drag : S/C cross section to drag [m^2]
%       .rc_cd   : S/C drag coefficient
%  - othr_params: Structure containing some physical parameters
%  - Env    : Structure containing the environment definition
% OUTPUTS:
%  - drag_pert: Vectorial acceleration due to drag in MEE2000 components 
%               [km/s^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------- START OF CODE -----------------------------
cbody_str  = Env.cbody;
cbody_spin = othr_params.cbody_spin;
Re         = othr_params.Re_cen;
% --------------------------- COMPUTE DENSITY -----------------------------
if ( strcmp(cbody_str,'EARTH') )
    if (Env.idrag == 1)
        % Scale height for the outer atmosphere density computation [km]
        H = 58.515;
        % Reference density at 400 km [kg/m^3]
        rho0 = 3.725*10^-12;
        alt0 = 400;
        % Altitude above a spherical Earth
        alt     = norm(state(1:3)) - Re;
        % Atmospheric density [kg/m^3]
        rho     = rho0 * exp( - (alt-alt0) / H);
    elseif (Env.idrag == 2)
        % Reference altitudes [km]
        alt0s = [150,180,200,250,300,350,400,450,500];
        % Reference densities [kg/m^3]
        rho0s = [2.07e-9, 5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, ...
                 9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13];
        % Scale heights [km]
        Hs    = [22.525,29.740,37.105,45.546,53.628,53.298,58.515,...
                 60.828,63.822];
        % Altitude above a spherical Earth
        alt     = norm(state(1:3)) - Re;
        ind     = 1;
        if (alt > alt0s(end) )
            ind = 10;
        elseif (alt < alt0s(1))
            ind = 2;
        else  
            while (alt > alt0s(ind) )
                ind = ind+1;
            end       
        end    
        % Define correct values for base altitude, density and scale height
        alt0 = alt0s(ind-1);
        rho0 = rho0s(ind-1);
        H    = Hs (ind-1);
        % Atmospheric density computation [kg/m^3]
        rho     = rho0 * exp( - (alt-alt0) / H);
    end
else
    rho = 0.0;
end

% Drag cross section
cs_drag = SC_params.cs_drag;
% Drag coefficient
cd      = SC_params.cd;

if (Env.idrag == 1)
    v_atm = 0;
elseif (Env.idrag == 2)  
    % Rotation matrix from MEE2000 to rotating reference frame
    cbody_rot_frame = ['IAU_',cbody_str];
    ROT_ROTFRAME_2_MEE2000 = cspice_pxform(cbody_rot_frame, 'J2000', et);
    omega_ver  = ROT_ROTFRAME_2_MEE2000 (:,3);
    omega      = omega_ver*cbody_spin;
    % Velocity of the atmosphere at the point of interest (in the inertial
    % frame)
    v_atm = cross (omega, state(1:3));
end

% S/C velocity relative to the atmosphere
Vrel      = state(4:6) - v_atm;
Vrel_mod  = norm(Vrel);
Vrel_vers = Vrel /Vrel_mod;

% Ballistic coefficient [kg/m^2]
BC        = state(7) / (cs_drag * cd);

% Acceleration due to drag [km/s^2]
drag_pert = - 1e3 * 0.5 * rho / BC * Vrel_mod^2 * Vrel_vers;


