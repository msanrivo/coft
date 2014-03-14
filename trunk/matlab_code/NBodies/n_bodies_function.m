%[This function integrates the ODE considering the effect of each body on the others.]
function dr = n_bodies_function(t,x,mu,n)

a=zeros(n*3,1);
for i=0:(n-1)
    for j=0:(n-1)
        if j~=i
        a(3*i+1)=a(3*i+1)-(mu(j+1)*(x(3*i+1)-x(3*j+1)))/(((x(3*i+1)-x(3*j+1)))^2+((x(3*i+2)-x(3*j+2)))^2+((x(3*i+3)-x(3*j+3)))^2)^(3/2);
        a(3*i+2)=a(3*i+2)-(mu(j+1)*(x(3*i+2)-x(3*j+2)))/(((x(3*i+1)-x(3*j+1)))^2+((x(3*i+2)-x(3*j+2)))^2+((x(3*i+3)-x(3*j+3)))^2)^(3/2);
        a(3*i+3)=a(3*i+3)-(mu(j+1)*(x(3*i+3)-x(3*j+3)))/(((x(3*i+1)-x(3*j+1)))^2+((x(3*i+2)-x(3*j+2)))^2+((x(3*i+3)-x(3*j+3)))^2)^(3/2);
        end
    end
end
        
dr=[x(3*n+1:end);a];

end