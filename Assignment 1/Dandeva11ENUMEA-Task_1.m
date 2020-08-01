clear
close all
clc

syms x y z v

z=(x^3 + (cos(y))/3)/(y-(sin(y))/2)

Tx=x/z*diff(z,x)

figure(1)
fsurf(Tx,[0,1,0,1])
xlabel('x')
ylabel('y')

Tx2= (3*x^3)/(x^3 + (cos(y))/3)
 
figure(2)
fsurf(Tx2,[0,1,0,1])
xlabel('x')
ylabel('y')


Ty=y/z*diff(z,y)

figure(3)
fsurf(Ty,[0,1,0,1])
xlabel('x')
ylabel('y')

Ty2=((-y/3)*sin(y))/(x^3+(cos(y))/3) - y/(y-0.5*sin(y)) + (0.5*y*cos(y))/(y-0.5*sin(y))
 
figure(4)
fsurf(Ty2,[0,1,0,1])
xlabel('x')
ylabel('y')
 

zs=subs(z,x^3,v)
Kpow=v/zs*diff(zs,v)
Kpow=subs(Kpow,v,x^3)
Kpow2=(x^3)/(x^3 + (cos(y))/3)

figure(5)
fsurf(Kpow,[0,1,0,1])
xlabel('x')
ylabel('y')

figure(6)
fsurf(Kpow2,[0,1,0,1])
xlabel('x')
ylabel('y')


zs=subs(z,cos(y),v)
Kcos=v/zs*diff(zs,v)
Kcos=subs(Kcos,v,cos(y))
Kcos2=((cos(y))/3)/(x^3 + (cos(y))/3)

figure(7)
fsurf(Kcos,[0,1,0,1])
xlabel('x')
ylabel('y')

figure(8)
fsurf(Kcos2,[0,1,0,1])
xlabel('x')
ylabel('y')


zs=subs(z,(1/3),v)
K3=v/zs*diff(zs,v)
K3=subs(K3,v,(1/3))
K3_2=((cos(y))/3)/(x^3 + (cos(y))/3)

figure(9)
fsurf(K3,[0,1,0,1])
xlabel('x')
ylabel('y')

figure(10)
fsurf(K3_2,[0,1,0,1])
xlabel('x')
ylabel('y')


zs=subs(z,(x^3 + (cos(y))/3),v)
Ksum=v/zs*diff(zs,v)
Ksum=subs(Ksum,v,(x^3 + (cos(y))/3))
Ksum2= 1

zs=subs(z,(y-(sin(y))/2),v)
Ksub=v/zs*diff(zs,v)
Ksub=subs(Ksub,v,(y-(sin(y))/2))
Ksub2= -1


zs=subs(z,(x^3 + (cos(y))/3)/(y-(sin(y))/2),v)
Kdiv=v/zs*diff(zs,v)
Kdiv=subs(Kdiv,v,(x^3 + (cos(y))/3)/(y-(sin(y))/2))
Kdiv2= 1

zs=subs(z,sin(y),v)
Ksin=v/zs*diff(zs,v)
Ksin=subs(Ksin,v,sin(y))
Ksin2=((sin(y))/2)/(y-(sin(y))/2)

figure(11)
fsurf(Ksin,[0,1,0,1])
xlabel('x')
ylabel('y')

figure(12)
fsurf(Ksin2,[0,1,0,1])
xlabel('x')
ylabel('y')


zs=subs(z,sin(y)/2,v*sin(y))
K2=v/zs*diff(zs,v)
K2=subs(K2,v,1/2)
K2_2=((sin(y))/2)/(y - (sin(y))/2)

figure(13)
fsurf(K2,[0,1,0,1])
xlabel('x')
ylabel('y')

figure(14)
fsurf(K2_2,[0,1,0,1])
xlabel('x')
ylabel('y')




