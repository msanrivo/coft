% Calcula la fecha y la hora a partir de los datos TLE
%
% function fecha=dia_hora(d,h) 
% 
% Author:   
% Date:     
 
% Entrada:
%       
%   d, dia juliano
%   h, fraccion de dia con 8 digitos
%	   
% Salida:
%
%	dia/mes
%	hh/mm/ss
%	segundos transcurridos en el dia


function tf=dia_hora(line_1);

year=floor((line_1(4))*0.001);

d=floor(line_1(4)-year*1000);

h=line_1(4)-year*1000-d;

if  d==0
    disp('31 diciembre')
elseif     d<=31
    dia=d;
    disp(dia)
    disp('enero')
elseif     d<=59
    dia=d-31;
    disp(dia)
    disp('febrero')
elseif    d<=90
    dia=d-59;
    disp(dia)
    disp('marzo')
elseif    d<=120
    dia=d-90;
    disp(dia)
    disp('abril')
elseif    d<=151
    dia=d-120;
    disp(dia)
    disp('mayo')
elseif    d<=181
    dia=d-151;
    disp(dia)
    disp('junio')
elseif    d<=212
    dia=d-181;
    disp(dia)
    disp('julio')
elseif    d<=243
    dia=d-212;
    disp(dia)
    disp('agosto')
elseif    d<=273
    dia=d-243;
    disp(dia)
    disp('septiembre')
elseif    d<=304
    dia=d-273;
    disp(dia)
    disp('octubre')
elseif    d<=334
    dia=d-304;
    disp(dia)
    disp('noviembre')
elseif    d<=364 
    dia=d-334;
    disp(dia)
    disp('diciembre')
end

d_1=d
a=h*24;
disp('hora');
disp(floor(a));
b=(a-floor(a))*60;
disp('minuto');
disp(floor(b));
c=(b-floor(b))*60;
disp('segundo');
disp(floor(c));

disp('un total de [segundos]')
disp(3600*floor(a)+60*floor(b)+floor(c))

fecha_1=(3600*floor(a)+60*floor(b)+floor(c));

line_3=input('line_3: EXAMPLE [1 25544 98067    08264.51782528  -.00002182  00000-0 -11606-4 0  2927] ');

year=floor((line_3(4))*0.001);

d=floor(line_3(4)-year*1000);

h=line_3(4)-year*1000-d;

if  d==0
    disp('31 diciembre')
elseif     d<=31
    dia=d;
    disp(dia)
    disp('enero')
elseif     d<=59
    dia=d-31;
    disp(dia)
    disp('febrero')
elseif    d<=90
    dia=d-59;
    disp(dia)
    disp('marzo')
elseif    d<=120
    dia=d-90;
    disp(dia)
    disp('abril')
elseif    d<=151
    dia=d-120;
    disp(dia)
    disp('mayo')
elseif    d<=181
    dia=d-151;
    disp(dia)
    disp('junio')
elseif    d<=212
    dia=d-181;
    disp(dia)
    disp('julio')
elseif    d<=243
    dia=d-212;
    disp(dia)
    disp('agosto')
elseif    d<=273
    dia=d-243;
    disp(dia)
    disp('septiembre')
elseif    d<=304
    dia=d-273;
    disp(dia)
    disp('octubre')
elseif    d<=334
    dia=d-304;
    disp(dia)
    disp('noviembre')
elseif    d<=364 
    dia=d-334;
    disp(dia)
    disp('diciembre')
end

d_2=d
a=h*24;
disp('hora');
disp(floor(a));
b=(a-floor(a))*60;
disp('minuto');
disp(floor(b));
c=(b-floor(b))*60;
disp('segundo');
disp(floor(c));

disp('un total de [segundos]')
disp(3600*floor(a)+60*floor(b)+floor(c))
fecha_2=(3600*floor(a)+60*floor(b)+floor(c));

tf=86400-fecha_1+fecha_2+(d_2-d_1-1)*86400;

end