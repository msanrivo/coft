%[This function adimensionalizes the problem dividing each magnitude by its characteristic magnitude.]

function [ M,mu,ro,vo,InitCond,tspan ] = adim(mc,Lc,tc,M,mu,ro,vo,tspan)

M=M/mc;
mu=mu/(Lc^3/tc^2);
ro=ro/Lc;
vo=vo/(Lc/tc);
InitCond=[ro,vo]';
tspan=tspan/tc;

end