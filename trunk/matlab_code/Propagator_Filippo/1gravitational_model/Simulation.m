function Simulation(t,X,R)

[xx,yy,zz] = sphere(100);
xx = xx.*R;
yy = yy.*R;
zz = zz.*R;

for i = 1 : length(t)
   
    surf(xx,yy,zz);
    hold on
    plot3(X(:,1),X(:,2),X(:,3));
    hold on
    plot3(X(i,1),X(i,2),X(i,3),'r+');
    hold off
    title(sprintf('Time passed: %1.1f', t(i)));
    
    drawnow;
    
    
end