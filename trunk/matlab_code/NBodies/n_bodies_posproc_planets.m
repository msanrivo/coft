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

%Coordinates of the baricenter: rg=sum(m*r)/sum(m) rg==wx3, w corresponds to the element, 3 to the coordinate
%Velocity of the baricenter: vg=sum(m*v)/sum(m) vg==wx3, w corresponds to the element, 3 to the coordinate
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

%Coordinates of each body w.r.t. the baricenter:
%rbg==wxnx3, w corresponds to the element, n to the body, and 3 to the coordinate
%Velocities of each body w.r.t. the baricenter:
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
s=[s;s];

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

%Earth-sat.:
a=4; %Central body
b=n; %Orbital body
figure(10),plot(T,normrr(:,b,a),'b',T,max(normrr(:,b,a)),'r',T,min(normrr(:,b,a)),'g'),title('Distance sat.-Earth')
xlabel('T')
ylabel('r')
figure(11),plot(T,normvr(:,b,a),'b',T,max(normvr(:,b,a)),'r',T,min(normvr(:,b,a)),'g'),title('Velocity sat.-Earth')
xlabel('T')
ylabel('v')