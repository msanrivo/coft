% Orbit radius crossing
function [lookfor stop direction] = orbit_radius_crossing(et,state,r_tg)

r =  sqrt(state(1)^2+state(2)^2+state(3)^2);

lookfor = r-r_tg;    %Searches for this expression set to 0
stop    = 1;         %Stop when event is located
direction = -1;       %Specifiy direction of motion at event