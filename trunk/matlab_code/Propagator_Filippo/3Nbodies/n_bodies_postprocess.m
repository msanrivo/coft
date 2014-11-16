%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script postprocesses the results. It plots graphs with the results
%
% Author: Pablo A. Machuca Varela
%
% Input:
% T              [length(T)x1]    Dimensionalized time
% X              [length(T)x6*n]  [x, y, z, vx, vy, vz] Dimensionalized results from the integration
% a              [1x1]            Central body of the two body problem for orbital parameters computation
% b              [1x1]            Orbital body of the two body problem for orbital parameters computation
%
% Output:
% rr             [length(T)xnxnx3] [xr, yr, zr]     Relative positions of each body w.r.t. the others. n corresponds to the body, n to the body with respect to whom it is been represented
% vr             [length(T)xnxnx3] [vxr, vyr, vzr]  Relative velocities of each body w.r.t. the others. n corresponds to the body, n to the body with respect to whom it is been represented
% normrr         [length(T)xnxn]   Relative distances of each body w.r.t. the others. n corresponds to the body, n to the body with respect to whom it is been represented
% normvr         [length(T)xnxn]   Norm of the relative velocities of each body w.r.t. the others. n corresponds to the body, n to the body with respect to whom it is been represented
% rg             [length(T)x3]     [xg, yg, zg]     Position of the barycenter 
% vg             [length(T)x3]     [vxg, vyg, vzg]  Velocity of the barycenter 
% normrg         [length(T)x1]     Norm of the barycenter position w.r.t. absolute reference system
% normvg         [length(T)x1]     Norm of the barycenter velocity w.r.t. absolute reference system
% rbg            [length(T)xnx3]   [xbg, ybg, zbg]      Positions w.r.t. the barycenter
% vbg            [length(T)xnx3]   [vxbg, vybg, vzbg]   Velocities w.r.t. the barycenter
% normrbg        [length(T)xn]     Relative distances w.r.t. the barycenter
% normvbg        [length(T)xn]     Norm of the relative velocities w.r.t. the barycenter
% an             [length(T)x1]     Semimajor axis of the orbit
% tn             [length(T)x1]     Theoretical period of the orbit
% meanan         [1x1]             Mean semimajor axis of the orbit
% meantn         [1x1]             Mean theoretical period of the orbit
% en             [length(T)x3]     [enx, eny, enz]    Eccentricity vector of the orbit
% costhn         [length(T)x1]     Cosine of the angle between the position vector and the initial position vector  
% normen         [length(T)x1]     Norm of the eccentricity vector of the orbit
% meannormen     [1x1]             Mean eccentricity of the orbit
% normrnp        [length(T)x1]     Periapsis of the orbit
% normrna        [length(T)x1]     Apoapsis of the orbit
% meannormrnp    [length(T)x1]     Mean periapsis of the orbit
% meannormrna    [length(T)x1]     Mean apoapsis of the orbit
% K              [length(T)x1]     Kinetic energy of the system
% U              [length(T)x1]     Potential energy of the system
% E              [length(T)x1]     Total energy of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Relative states of each body w.r.t. the others
rr = zeros([length(T),n,n,3]);
vr = zeros([length(T),n,n,3]);
for i=1:n
    for j=1:n
        rr(:,i,j,1) = X(:,3*(i-1)+1)-X(:,3*(j-1)+1);
        rr(:,i,j,2) = X(:,3*(i-1)+2)-X(:,3*(j-1)+2);
        rr(:,i,j,3) = X(:,3*(i-1)+3)-X(:,3*(j-1)+3);
        vr(:,i,j,1) = X(:,3*n+3*(i-1)+1)-X(:,3*n+3*(j-1)+1);
        vr(:,i,j,2) = X(:,3*n+3*(i-1)+2)-X(:,3*n+3*(j-1)+2);
        vr(:,i,j,3) = X(:,3*n+3*(i-1)+3)-X(:,3*n+3*(j-1)+3);
    end
end
%
normrr = sqrt(sum(rr.^2,4));
normvr = sqrt(sum(vr.^2,4));
%
%State of the barycenter
rg = zeros([length(T),3]);
vg = zeros([length(T),3]);
for i=0:n-1
    rg(:,1) = rg(:,1) + M(i+1)*X(:,3*i+1)/sum(M);     %rg = sum(m*r)/sum(m)
    rg(:,2) = rg(:,2) + M(i+1)*X(:,3*i+2)/sum(M);
    rg(:,3) = rg(:,3) + M(i+1)*X(:,3*i+3)/sum(M);
    vg(:,1) = vg(:,1) + M(i+1)*X(:,3*n+3*i+1)/sum(M); %vg = sum(m*v)/sum(m)
    vg(:,2) = vg(:,2) + M(i+1)*X(:,3*n+3*i+2)/sum(M);
    vg(:,3) = vg(:,3) + M(i+1)*X(:,3*n+3*i+3)/sum(M);
end
rg(:,1) = rg(:,1)-rg(1,1);
rg(:,2) = rg(:,2)-rg(1,2);
rg(:,3) = rg(:,3)-rg(1,3);
%
normrg = sqrt(sum(rg.^2,2));
normvg = sqrt(sum(vg.^2,2));
%
%Relative states w.r.t. the barycenter
rbg = zeros([length(T),n,3]);
vbg = zeros([length(T),n,3]);
for i=0:n-1
    rbg(:,i+1,1) = X(:,3*i+1)-rg(:,1);
    rbg(:,i+1,2) = X(:,3*i+2)-rg(:,2);
    rbg(:,i+1,3) = X(:,3*i+3)-rg(:,3);
    vbg(:,i+1,1) = X(:,3*n+3*i+1)-vg(:,1);
    vbg(:,i+1,2) = X(:,3*n+3*i+2)-vg(:,2);
    vbg(:,i+1,3) = X(:,3*n+3*i+3)-vg(:,3);
end
%
normrbg = sqrt(sum(rbg.^2,3));
normvbg = sqrt(sum(vbg.^2,3));
%
%Orbital parameters
%
%%%%%%%%%%%
%Sun-Earth
%Central body and orbital body
a = 1;
b = 4; 
%
%Semimajor axis and period
a3     = (-normvr(:,b,a).^2/mu(a)+2./normrr(:,b,a)).^(-1); %Energy equation: 1/2*v^2-mu/r=-mu/(2*a)
t3    =  2*pi*a3.^(3/2)/((mu(a))^(1/2));                   %Third Kepler's Law: T^2 = 4*pi^2/(G*M)*a^3
meana3 = mean(a3);
meant3 = mean(t3);
%
%Eccentricity and cosine
e3       = zeros([length(T),3]);
costh3   = zeros([length(T),1]); 
for i=1:length(T)
    e3(i,:)   = ((normvr(i,b,a)^2-mu(a)/normrr(i,b,a))*rr(i,b,a,:)-dot(rr(i,b,a,:),vr(i,b,a,:))*vr(i,b,a,:))/mu(a); %e=((v^2-mu/r)*r - (r·v)*v)/mu
    costh3(i) = dot(rr(i,b,a,:),rr(1,b,a,:))./(normrr(i,b,a)*normrr(1,b,a));                                        %r·r0=|r||r0|*cos(theta)
end
%
norme3       = sqrt(sum(e3.^2,2));
meannorme3   = mean(norme3);
%
%Periapsis and apoapsis
normr3p       = a3.*(1-norme3); %rp = a*(1-e) 
normr3a       = a3.*(1+norme3); %ra = a*(1+e)
meannormr3p   = mean(normr3p);
meannormr3a   = mean(normr3a);
%
%%%%%%%%%%%
%sat.-Earth
%Central body and orbital body
a = 4;
b = 5;
%
%Semimajor axis and period
a313  = (-normvr(:,b,a).^2/mu(a)+2./normrr(:,b,a)).^(-1); %Energy equation: 1/2*v^2-mu/r=-mu/(2*a)
t313 = 2*pi*a313.^(3/2)/((mu(a))^(1/2));                  %Third Kepler's Law: T^2 = 4*pi^2/(G*M)*r^3
meana313  = mean(a313);
meant313 = mean(t313);
%
%Eccentricity and cosine
e313     = zeros([length(T),3]);
costh313 = zeros([length(T),1]); 
for i=1:length(T)
    e313(i,:)   = ((normvr(i,b,a)^2-mu(a)/normrr(i,b,a))*rr(i,b,a,:)-dot(rr(i,b,a,:),vr(i,b,a,:))*vr(i,b,a,:))/mu(a); %e=((v^2-mu/r)*r - (r·v)*v)/mu
    costh313(i) = dot(rr(i,b,a,:),rr(1,b,a,:))./(normrr(i,b,a)*normrr(1,b,a));                                        %r·r0=|r||r0|*cos(theta)
end
%
norme313     = sqrt(sum(e313.^2,2)); %Eccentricity sat-Earth
meannorme313 = mean(norme313);
%
%Periapsis and apoapsis
normr313p = a313.*(1-norme313); %rp = a*(1-e) 
normr313a = a313.*(1+norme313); %ra = a*(1+e)
meannormr313p = mean(normr313p);
meannormr313a = mean(normr313a);
%
%%%%%%%%%%%
%Energy of the system 
%
%Kinetic energy
K = zeros([length(T),1]);
for i=1:n
    K = K + (1/2)*M(i)*normvbg(:,i).^2;        %K = 1/2*m*v^2
end
%
%Potential energy
U = zeros([length(T),1]);
for i=1:(n-1)
    for j=i:(n-1)
        U = U - mu(i)*M(j+1)./normrr(:,i,j+1); %U = -mu*M/|r2-r1|
    end
end
%
E = K + U;
%
%Plots
%
%To plot every body in a different color
s = ['b'; 'g'; 'r'; 'c'; 'm'; 'y'; 'k'];
s = [s; s; s];
for i=0:(n-1)
    figure(1)
    plot3(X(:,3*i+1),X(:,3*i+2),X(:,3*i+3),s(i+1)) %3-D plot
    xlabel('x(m)'), ylabel('y(m)'), zlabel('z(m)')
    hold on
end
hold off
%
figure(2)
plot(T,E,'b') %Energy of the system
title('Energy of the system')
xlabel('T(s)'), ylabel('E(J)')
%
%Earth-Sun:
figure(3)
plot(T,costh3) %Cosine
title(['cos(theta) ' num2str(a) '-' num2str(b)])
xlabel('T'), ylabel('cos(theta)')
%
%Central body (a) and orbital body (b)
a=1; 
b=4; 
%
figure(4)
plot(T,normrr(:,b,a),'b',T,max(normrr(:,b,a)),'r',T,min(normrr(:,b,a)),'g') %Distance a-b
title(['Distance ' num2str(a) '-' num2str(b)])
xlabel('T'), ylabel('r')
%
figure(5)
plot(T,normvr(:,b,a),'b',T,max(normvr(:,b,a)),'r',T,min(normvr(:,b,a)),'g') %Velocity magnitude a-b
title(['Velocity ' num2str(a) '-' num2str(b)])
xlabel('T'), ylabel('v')
