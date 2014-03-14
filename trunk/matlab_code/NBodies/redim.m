%[This function redimensionalizes the problem multiplying each magnitude by its characteristic magnitude.]

function [ M,mu,X,T ] = redim(mc,Lc,tc,M,mu,X,T,n)

M=M*mc;
mu=mu*(Lc^3/tc^2);
X(:,1:3*n)=X(:,1:3*n)*Lc;
X(:,(3*n+1):end)=X(:,(3*n+1):end)*Lc/tc;
T=T*tc;

end