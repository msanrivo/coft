function [lookfor stop direction] = node_crossing_detection(et,state)

z_coordinate =  state(3);

lookfor = z_coordinate ; %Searches for this expression set to 0
stop    = 1;     %Stop when event is located
direction = 1;  %Specifiy direction of motion at event