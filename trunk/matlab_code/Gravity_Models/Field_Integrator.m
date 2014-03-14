function state_output = Field_Integrator(t, X, GM,R,itggoce02_cell,NN,tf)

 cart_position = [X(1);X(2);X(3)];
 v = [X(4);X(5);X(6)];
 
 [lambda,phi,r] = cart2sph(X(1),X(2),X(3));
 
 rdotdot = Acceleration_Field3(r,phi,lambda,GM,R,itggoce02_cell,NN);
 d2r_fict = NonInertialTerms(v,cart_position);
 
 state_output = [v;(rdotdot-d2r_fict)];
    
 clc
 fprintf('Progress %1.2f percent \n',t/tf*100);
    
end

