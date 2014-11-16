function  d2r_fict = FictAcc(v,r)

%{

    Non-Inertial Acceleration In Earth Fixed Frame due to Earth Spin

    Author: Diego García Pardo (UNIVERSITY CARLOS III OF MADRID)

    As the acceleration must be presented in cartesian coordinates, 
    no especial considerations about gradients must be taken into account

    It is assumed that the Earth Rotation Speed is Constant
    It is assumed that the Earth Rotation Speed is w = 7.2921150 e - 5 rads/s 

    v & r must be column vectors
    
%}

w = 7.2921150e-5;

coriolis = 2.*[-v(2)*w ; v(1)*w ; 0];
centrifugal = [-r(1)*w^2; -r(2)*w^2; 0];

d2r_fict = coriolis + centrifugal;


end