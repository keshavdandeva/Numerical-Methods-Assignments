
clear 
close all
clc


d = de2bi(0:1023);
d(d==0) = -1;
d = 10^(-12)*d;

vec = zeros(1,10);

 for i = 1:1024
     for j = 1:10
         vec(j) = d(i,j);
     end
              
     z =@(x,y) abs(((x^3.*(1+ vec(1))^3.*(1+vec(2)) + (1./3).*cos(y.*(1+vec(3))).*(1+ vec(4)).*(1+ vec(5))).*(1+ vec(6)).*(1 + vec(7)))./((y.*(1+ vec(3)) - (1./2).*sin(y.*(1+ vec(3))).*(1+vec(8)).*(1+vec(9))).*(1+vec(10))));
     
     x = linspace(1,10,500);
     y = linspace(1,10,500);

     [X,Y] = meshgrid(x,y);
     sh = surf(X, Y, z(X,Y));

     zd=get(sh,'zdata');
     zmax=max(max(zd));
     m(i) = zmax;
     zs =@(x,y) (x^3 + (cos(y))./3)./(y-(sin(y))./2);
     x = linspace(1,10,500);
     y = linspace(1,10,500);

     [X,Y] = meshgrid(x,y);
     sh = surf(X, Y, zs(X,Y));

     zd=get(sh,'zdata');
     zmax=max(max(zd));
     expected_z = zmax;

     relative_error(i) = abs((m(i) - expected_z)./expected_z);
 end

max(relative_error)


 