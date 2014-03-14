function [vx, vy, vz] = cart_vel(r,phi,th,drdt, dphidt, dthdt)

vx = drdt * cos(phi)*cos(th) - r*sin(phi)*cos(th)*dphidt - r*cos(phi)*sin(th)*dthdt;
vy = drdt*cos(phi)*sin(th) - r*sin(phi)*sin(th)*dphidt  + r*cos(phi)*cos(th)*dthdt;
vz = drdt*sin(phi) + r*cos(phi)*dphidt;


end

