% Simple solver for Kepler's equation
% M = E - e*sin(E)
% MMM20141001


M = 235.4/180*pi;
e = 0.4;
E = M;
max_iter = 200;
tol = 1e-8;

for i=0:max_iter
    err = E - e*sin(E) - M;    
    deltaE = err/(1-e*cos(E));
    fprintf('Iter %i: E = %e deg, error = %e, deltaE/E = %e\n',i,E*180/pi,err,deltaE/E)    
    if abs(err) < tol && abs(deltaE/E) < tol
        fprintf('Convergence reached on iter %i: E = %e deg, error = %e\n',i,E,err)
        break
    end
    E = E - deltaE;
end