function state_output = Field_Integrator(t, X, GM,R,GeoP_Coeff,NN,tf)

%{

    Kinematic & Dynamic Relationships To Produce Vector Of Velocities And
    Accelerations

    Author: Diego García Pardo (UNIVERSITY CARLOS III OF MADRID)

    Input Definitions:

        -- t            : Time Instant (Integrator Time)                                    [s]
        -- X            : Vector Of Positions And Velocities At Previous Time Step          [m;m/s]         
        -- GM           : Standard Gravitational Constant                                   [m^3 s^-2]
        -- GeoP_Coeff   : Geopotential Coefficientes For the Associated
                          Planet (FULLY NORMALIZED)                                         
        -- NN           : Model Size                                                        []
        -- tf           : Final Time Of integration (use for purposes of progress)          [s]


%}

%% Computing Position And Velocity At Given Instant

cart_position = [X(1);X(2);X(3)];
v = [X(4);X(5);X(6)];

%% Transforming Cartesian Coordinates to Spherical Coordinates (Required for Acceleartion Field Algorithm)

[lambda,phi,r] = cart2sph(X(1),X(2),X(3));

%% Computing Acceleration Field And Non Inertial Terms

rdotdot = Acc_Field(r,phi,lambda,GM,R,GeoP_Coeff,NN);  % Acceleration
% d2r_fict = FictAcc(v,cart_position);                   % Non Inertial Terms
d2r_fict =0;
%% Outputting Results

state_output = [v;(rdotdot-d2r_fict)];                 % Column Vector Of velocities And Accelerations

clc
fprintf('Progress %1.2f percent \n',t/tf*100);         % Progress Percentage (CAUTION: MAY BE SOURCE OF SLOWNESS)

end

