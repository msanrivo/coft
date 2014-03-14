%[This script defines the constant G, data of the bodies, initial conditions and tspan.]
%[It calls to the adimensionalization, integration, redimensionalization and posprocessing functions.]

clear all
close all
clc
%%%%%%%%%%%%
octave = 0;
%%%%%%%%%%%% 1 AU=149597870691m 1d=86400s
%Data:
G=1.4878e-034; %1.4878e-034 [AU^3/(d^2*kg)] 6.67259*10^-11 [m^3/s^2*kg]
%0==Sun 1==Mercury 2==Venus 3==Earth 4==Mars 5==Jupiter 6==Saturn 7==Uranus 8==Neptune 9==Pluto 31==Moon
m0=1.9891*10^30; %Mass of the Sun
m1=3.302e+023;
m2=4.8685e+024;
m3=5.97219*10^24; %kg Mass of the Earth
m31=7.349*10^22; %kg Mass of the satellite
m4=6.4191e+023;
m5=1.8992e+027; %kg Mass of Jupiter
m6=5.6866e+026;
m7=8.6849e+025;
m8=1.0247e+026;
m9=1.4734e+022;
mDA=40000000;
M=[m0;m1;m2;m3;m31;m4;m5;m6;m7;m8;m9;mDA];
mu=G*M;
n=numel(M);
%%%%%%%%%%%%
%Initial conditions: 2013-02-15 Origin: "Barycenter" Reference plane: "Ecliptic and mean equinox of reference epoch"
%Sun:
r0o=[-1.006136334109344E-03 -2.504882138940138E-03 -4.804565344486537E-05]; %AU
v0o=[6.222559730391206E-06 -8.904765292268053E-07 -1.354277838835098E-07]; %AU/d
%Mercury:
r1o=[1.340344577471494E-01  2.751308390697875E-01  1.024619965349763E-02]; %AU
v1o=[-3.091640845373111E-02  1.340631715635793E-02  3.932533766809858E-03]; %AU/d
%Venus:
r2o=[3.832056977407786E-01 -6.202928528888464E-01 -3.068724150838091E-02]; %AU
v2o=[1.704643655404703E-02  1.061031081630652E-02 -8.381812228809199E-04]; %AU/d
%Earth:
r3o=[-8.226211107254137E-01  5.457906263566129E-01 -6.519569880195386E-05]; %AU
v3o=[-9.821543261542958E-03 -1.438149210007801E-02  8.498810918980957E-07]; %AU/d
%Moon:
r31o=[-8.202317465636982E-01  5.468822579815011E-01  3.611907788126580E-05]; %AU
v31o=[-1.002905393392062E-02 -1.384582117054919E-02 -4.335204885740597E-05]; %AU/d
%Mars:
r4o=[1.361357874098748E+00 -2.480646841875632E-01 -3.864253863458948E-02]; %AU
v4o=[3.024812856420692E-03  1.496466333928509E-02  2.393074753960224E-04]; %AU/d
%Jupiter:
r5o=[1.090982664998702E+00  4.956614946535135E+00 -4.507977442301758E-02]; %AU
v5o=[-7.461377241667266E-03  1.982418921945255E-03  1.587313253301721E-04]; %AU/d
%Saturn:
r6o=[-7.951735276643929E+00 -5.723294868705614E+00  4.159710773510763E-01]; %AU
v6o=[2.956251284722622E-03 -4.541265150103193E-03 -3.867530584571000E-05]; %AU/d
%Uranus:
r7o=[1.986432265405680E+01  2.738096414520781E+00 -2.471803929877430E-01]; %AU
v7o=[-5.658304707907710E-04  3.712922858232250E-03  2.112521676132228E-05]; %AU/d
%Neptune:
r8o=[2.662393385906048E+01 -1.379954789611723E+01 -3.293920251777657E-01]; %AU
v8o=[1.423699074369510E-03  2.805614721239486E-03 -9.058685879265006E-05]; %AU/d
%Pluto:
r9o=[5.249904789126634E+00 -3.190348331761015E+01  1.895275741923878E+00]; %AU
v9o=[3.157333207199485E-03 -1.180187735535159E-04 -9.006600044110265E-04]; %AU/d

rDAo=[-8.218648861333863E-01  5.447437160747355E-01 -2.806397970305048E-03];
vDAo=[-1.105039092255858E-02 -1.313424560990138E-02  3.179724042324179E-03];

ro=[r0o,r1o,r2o,r3o,r31o,r4o,r5o,r6o,r7o,r8o,r9o,rDAo];
vo=[v0o,v1o,v2o,v3o,v31o,v4o,v5o,v6o,v7o,v8o,v9o,vDAo];
%%%%%%%%%%%%%%
tspan=[0 (2*365)];
%tspan=[0:1:365*2];
%%%%%%%%%%%%%%
%Adimensionalization:
b=7; %s=Chosen body for characteristic mass, length and time
mc=M(b); %Characteristic mass
Lc=norm(ro((3*(b-1)+1):(3*(b-1)+3))); %Characteristic length. Initial distance Earth-barycenter
tc=2*pi*norm(ro((3*(b-1)+1):(3*(b-1)+3)))^(3/2)/((mu(1))^(1/2)); %Characteristic time. Theoretical period of the Earth. Kepler's 3rd law

[ M,mu,ro,vo,InitCond,tspan ] = adim(mc,Lc,tc,M,mu,ro,vo,tspan);
%%%%%%%%%%%%%%
%Integration:
% Octave
if (octave)
tic
tspan_Oc = linspace(tspan(1), tspan(end)); 
X=lsode(@(x,t)n_bodies_function_Oc(x,t,mu,n), InitCond, tspan_Oc);
toc
T=tspan_Oc; 
else
%%%%%%%%%%%%%%
% Matlab
options=odeset('AbsTol',1e-14,'Reltol',1e-11);
tic
[T X]=ode45(@(t,x)n_bodies_function(t,x,mu,n), tspan, InitCond, options);
toc
end
%%%%%%%%%%%%%%
%Results:
%%%%%%%%%%%%%%
%Redimensionalization: %1 AU=149597870691m 1d=86400s
mc=mc*1; %Characteristic time multiplied by the conversion factor
Lc=Lc*149597870691; %Characteristic length multiplied by the conversion factor
tc=tc*86400; %Characteristic time multiplied by the conversion factor

[ M,mu,X,T ] = redim(mc,Lc,tc,M,mu,X,T,n);
%%%%%%%%%%%%%%
%Posprocessing:
n_bodies_posproc;
%End.