clear
close all
clc

syms x y z

z =@(x,y) (abs((3.*x^3)./(x^3 + (cos(y))./3)) + abs(((-y./3).*sin(y))./(x^3+(cos(y))./3) - y./(y-(0.5).*sin(y)) + ((0.5).*y.*cos(y))./(y-(0.5).*sin(y))) + abs((x^3)./(x^3 + (cos(y))./3)) + abs(((cos(y))./3)./(x^3 + (cos(y))./3)) + abs(((cos(y))./3)./(x^3 + (cos(y))./3)) + 3 + abs(((sin(y))./2)./(y-(sin(y))./2)) + abs(((sin(y))./2)./(y - (sin(y))./2)));

x = linspace(1,10,1000);
y = linspace(1,10,1000);

[X,Y] = meshgrid(x,y);
sh = surf(X, Y, z(X,Y));

zd=get(sh,'zdata');
zmax=max(max(zd))