function value = legendre(P,l,m,gamma)

lat = asin(gamma);
ll = l - 1;
mm = m - 1;

if (mm == 0)
    value = ( (2*ll-1)*gamma*P(l-1,1)-(ll-1)*P(l-2,1) ) / ll;
elseif (mm >0 && mm < ll)
    value = P(l-2,m) + (2*ll-1) * cos(lat) * P(l-1,m-1);
elseif (mm == ll)
    value = (2*ll-1) * cos(lat) * P(l-1,l-1);
elseif (mm > ll)
    value = 0;
end