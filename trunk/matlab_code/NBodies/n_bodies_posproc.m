%[This script posprocesses the data after the integration.]
%[Expresses each body w.r.t. the others, w.r.t. the barycenter,
%calculates the energy of the system and calculates some orbital parameters
%such as the semimajor axis, periapsis, apoapsis and eccentricity vector.]
%[It also plots the results.]

%X==[x1 y1 z1 x2 y2 z2 ... vx1 vy1 vz1 vx2 vy2 vz2 ...]
w=numel(T);
%Relative states of each body w.r.t. the others
%rr==matrix wxnxnx3. w corresponds to the element, n to the body, n to the other body with respect to whom it is been represented, and 3 to the coordinate.
%vr==matrix wxnxnx3. w corresponds to the element, n to the body, n to the other body with respect to whom it is been represented, and 3 to the coordinate.
rr=zeros(w,n,n,3);
vr=zeros(w,n,n,3);
for i=1:n
    for j=1:n
        rr(:,i,j,1)=X(:,3*(i-1)+1)-X(:,3*(j-1)+1);
        rr(:,i,j,2)=X(:,3*(i-1)+2)-X(:,3*(j-1)+2);
        rr(:,i,j,3)=X(:,3*(i-1)+3)-X(:,3*(j-1)+3);
        vr(:,i,j,1)=X(:,3*n+3*(i-1)+1)-X(:,3*n+3*(j-1)+1);
        vr(:,i,j,2)=X(:,3*n+3*(i-1)+2)-X(:,3*n+3*(j-1)+2);
        vr(:,i,j,3)=X(:,3*n+3*(i-1)+3)-X(:,3*n+3*(j-1)+3);
    end
end
normrr=sqrt(sum(rr.^2,4)); %Relative distance: wxnxn
normvr=sqrt(sum(vr.^2,4)); %Relative velocity: wxnxn

%Coordinates of the barycenter: rg=sum(m*r)/sum(m) rg==wx3, w corresponds to the element, 3 to the coordinate
%Velocity of the barycenter: vg=sum(m*v)/sum(m) vg==wx3, w corresponds to the element, 3 to the coordinate
rg=zeros(w,3);
vg=zeros(w,3);
for i=0:n-1
    rg(:,1)=rg(:,1)+M(i+1)*X(:,3*i+1)/sum(M);
    rg(:,2)=rg(:,2)+M(i+1)*X(:,3*i+2)/sum(M);
    rg(:,3)=rg(:,3)+M(i+1)*X(:,3*i+3)/sum(M);
    vg(:,1)=vg(:,1)+M(i+1)*X(:,3*n+3*i+1)/sum(M);
    vg(:,2)=vg(:,2)+M(i+1)*X(:,3*n+3*i+2)/sum(M);
    vg(:,3)=vg(:,3)+M(i+1)*X(:,3*n+3*i+3)/sum(M);
end
rg(:,1)=rg(:,1)-rg(1,1);
rg(:,2)=rg(:,2)-rg(1,2);
rg(:,3)=rg(:,3)-rg(1,3);
normvg=sqrt(vg(:,1).^2+vg(:,2).^2+vg(:,3).^2); %Modulus of vg

%Coordinates of each body w.r.t. the barycenter:
%rbg==wxnx3, w corresponds to the element, n to the body, and 3 to the coordinate
%Velocities of each body w.r.t. the barycenter:
%vbg==wxnx3, w corresponds to the element, n to the body, and 3 to the coordinate
rbg=zeros(w,n,3);
vbg=zeros(w,n,3);
for i=0:n-1
rbg(:,i+1,1)=X(:,3*i+1)-rg(:,1);
rbg(:,i+1,2)=X(:,3*i+2)-rg(:,2);
rbg(:,i+1,3)=X(:,3*i+3)-rg(:,3);
vbg(:,i+1,1)=X(:,3*n+3*i+1)-vg(:,1);
vbg(:,i+1,2)=X(:,3*n+3*i+2)-vg(:,2);
vbg(:,i+1,3)=X(:,3*n+3*i+3)-vg(:,3);
end
normrbg=sqrt(sum(rbg.^2,3)); %Relative distance: wxnxn
normvbg=sqrt(sum(vbg.^2,3)); %Relative velocity: wxnxn

%Orbital parameters:
%Semimajor axis. Energy equation: 1/2*v^2-mu/r=-mu/(2*a)
%Period. Kepler's 3rd law
%Sun-Earth
a=1; %Central body
b=4; %Orbital body
a3=(-normvr(:,b,a).^2/mu(a)+2./normrr(:,b,a)).^(-1);
tc3=2*pi*a3.^(3/2)/((mu(a))^(1/2));
%sat.-Earth
a=4; %Central body
b=5; %Orbital body
a313=(-normvr(:,b,a).^2/mu(a)+2./normrr(:,b,a)).^(-1);
tc313=2*pi*a313.^(3/2)/((mu(a))^(1/2));
%Sun-Jupiter
a=1; %Central body
b=7; %Orbital body
a5=(-normvr(:,b,a).^2/mu(a)+2./normrr(:,b,a)).^(-1);
tc5=2*pi*a5.^(3/2)/((mu(a))^(1/2));

%Cosine of the angle between the position vector and the initial position vector: r·r0=|r||r0|*cos(theta)
%Eccentricities: e=((v^2-mu/r)*r - (r·v)*v)/mu
costh3=zeros(w,1);
costh313=zeros(w,1);
costh5=zeros(w,1);
e3=zeros(w,3);
e313=zeros(w,3);
e5=zeros(w,3);
for i=1:w
%Earth-Sun orbit
a=1; %Central body
b=4; %Orbital body
costh3(i)=dot(rr(i,b,a,:),rr(1,b,a,:))./(normrr(i,b,a)*normrr(1,b,a));
e3(i,:)=((normvr(i,b,a)^2-mu(a)/normrr(i,b,a))*rr(i,b,a,:)-dot(rr(i,b,a,:),vr(i,b,a,:))*vr(i,b,a,:))/mu(a);
%Moon-Earth orbit
a=4; %Central body
b=5; %Orbital body
costh313(i)=dot(rr(i,b,a,:),rr(1,b,a,:))./(normrr(i,b,a)*normrr(1,b,a));
e313(i,:)=((normvr(i,b,a)^2-mu(a)/normrr(i,b,a))*rr(i,b,a,:)-dot(rr(i,b,a,:),vr(i,b,a,:))*vr(i,b,a,:))/mu(a);
%Jupiter-Sun orbit
a=1; %Central body
b=7; %Orbital body
costh5(i)=dot(rr(i,b,a,:),rr(1,b,a,:))./(normrr(i,b,a)*normrr(1,b,a));
e5(i,:)=((normvr(i,b,a)^2-mu(a)/normrr(i,b,a))*rr(i,b,a,:)-dot(rr(i,b,a,:),vr(i,b,a,:))*vr(i,b,a,:))/mu(a);
end
norme3=sqrt(sum(e3.^2,2)); %Eccentricity Sun-Earth
norme313=sqrt(sum(e313.^2,2)); %Eccentricity sat-Earth
norme5=sqrt(sum(e5.^2,2)); %Eccentricity Sun-Jup

%Periapsis and apoapsis: rp=a*(1-e) ra=a*(1+e)
normr3p=a3.*(1-norme3); %Perihelion Earth
normr3a=a3.*(1+norme3); %Aphelion Earth
normr313p=a313.*(1-norme313); %Perihelion sat
normr313a=a313.*(1+norme313); %Aphelion sat
normr5p=a5.*(1-norme5); %Perihelion Jup
normr5a=a5.*(1+norme5); %Aphelion Jup

%Energy of the system: 
K=zeros(w,1);
for i=1:n
    K=K + (1/2)*M(i)*normvbg(:,i).^2; %Kinetic energy
end
U=zeros(w,1);
for i=1:(n-1)
    for j=i:(n-1)
        U=U - mu(i)*M(j+1)./normrr(:,i,j+1); %Potential energy
    end
end
E=K + U;

%To plot every body in a different color
s=['b';'g';'r';'c';'m';'y';'k'];
s=[s;s;s];

for i=0:(n-1)
figure(1),plot3(X(:,3*i+1),X(:,3*i+2),X(:,3*i+3),s(i+1)) %3-D plot
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
hold on
end
hold off
figure(2),plot(T,E,'b'),title('Energy of the system')
xlabel('T(s)')
ylabel('E(J)')

%Earth-Sun:
figure(3),plot(T,costh3),title('cos(theta) Sun-Earth')
xlabel('T')
ylabel('cos(theta)')
a=1; %Central body
b=4; %Orbital body
figure(4),plot(T,normrr(:,b,a),'b',T,max(normrr(:,b,a)),'r',T,min(normrr(:,b,a)),'g'),title('Distance Sun-Earth')
xlabel('T')
ylabel('r')
figure(5),plot(T,normvr(:,b,a),'b',T,max(normvr(:,b,a)),'r',T,min(normvr(:,b,a)),'g'),title('Velocity Sun-Earth')
xlabel('T')
ylabel('v')

%Jup-Sun:
figure(6),plot(T,costh5),title('cos(theta) Sun-Jup')
xlabel('T')
ylabel('cos(theta)')
a=1; %Central body
b=7; %Orbital body
figure(7),plot(T,normrr(:,b,a),'b',T,max(normrr(:,b,a)),'r',T,min(normrr(:,b,a)),'g'),title('Distance Sun-Jup')
xlabel('T')
ylabel('r')
figure(8),plot(T,normvr(:,b,a),'b',T,max(normvr(:,b,a)),'r',T,min(normvr(:,b,a)),'g'),title('Velocity Sun-Jup')
xlabel('T')
ylabel('v')

%Earth-sat.:
figure(9),plot(T,costh313),title('cos(theta) sat-Earth')
xlabel('T')
ylabel('cos(theta)')
a=4; %Central body
b=5; %Orbital body
figure(10),plot(T,normrr(:,b,a),'b',T,max(normrr(:,b,a)),'r',T,min(normrr(:,b,a)),'g'),title('Distance sat.-Earth')
xlabel('T')
ylabel('r')
figure(11),plot(T,normvr(:,b,a),'b',T,max(normvr(:,b,a)),'r',T,min(normvr(:,b,a)),'g'),title('Velocity sat.-Earth')
xlabel('T')
ylabel('v')

disp(['Period of the orbit of the Earth= ' num2str(mean(tc3)) 's = ' num2str(mean(tc3)/(3600*24)) ' days.'])
disp(['Semimajor axis of the orbit of the Earth= ' num2str(mean(a3)) 'm.'])
disp(['Eccentricity of the orbit of the Earth= ' num2str(mean(norme3)) '.'])
disp(['Periapsis E.= ' num2str(mean(normr3p)) 'm.'])
disp(['Apoapsis E= ' num2str(mean(normr3a)) 'm.'])

disp(['Period of the orbit of the sat.= ' num2str(mean(tc313)) 's = ' num2str(mean(tc313)/(3600*24)) ' years.'])
disp(['Semimajor axis of the orbit sat.= ' num2str(mean(a313)) 'm.'])
disp(['Eccentricity of the orbit sat.= ' num2str(mean(norme313)) '.'])
disp(['Periapsis sat.= ' num2str(mean(normr313p)) 'm.'])
disp(['Apoapsis sat= ' num2str(mean(normr313a)) 'm.'])

disp(['Period of the orbit of Jupiter= ' num2str(mean(tc5)) 's = ' num2str(mean(tc5)/(3600*24)) ' days.'])
disp(['Semimajor axis of the orbit of Jupiter= ' num2str(mean(a5)) 'm.'])
disp(['Eccentricity of the orbit of Jupiter= ' num2str(mean(norme5)) '.'])
disp(['Periapsis Jupiter= ' num2str(mean(normr5p)) 'm.'])
disp(['Apoapsis Jupiter= ' num2str(mean(normr5a)) 'm.'])