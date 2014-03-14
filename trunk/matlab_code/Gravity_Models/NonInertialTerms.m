function  d2r_fict = NonInertialTerms(v,r)

%{

    In this function, we output the value of the non-inertial accelerations
    As the earth moves around itself

    As the acceleration must be presented in cartesian coordinates, 
    no especial considerations about gradients must be taken into account

    It is assumed that the Earth Rotation Speed is Constant
    It is assumed that the Earth Rotation Speed is w = 7.2921150 e - 5 rads/s 

    v & r must be column vectors
    
%}

w =[0;0;7.2921150e-5];

coriolis = cross(w,v);
centrifugal = cross(w,cross(w,r));

d2r_fict = coriolis + centrifugal;


end