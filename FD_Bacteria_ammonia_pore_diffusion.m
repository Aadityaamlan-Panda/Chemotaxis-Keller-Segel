clear all
close all
clc
L=1;
T=1;
N=201; %space
M=201; %time
% Db=linspace(10^-13,10^-12,100);
Db=10^-13;
Dc=10^-9;
X0=10^-11;
% X0=linspace(10^-11,10^-9,100);
dx=L/(N-1);
x=linspace(0,1,N);
t=linspace(0,1,M);
dt=T/(M-1);
b=zeros(N,M);
b(:,1)=1;
l=[pi/2 3*pi/2 5*pi/2];
A=4.*(cos(l)-1)./(2.*l-sin(2.*l));

for n=1:M-1
    for i=2:N-1
       c1 = (1+exp(-l(1)^2*t(n))*A(1)*sin(l(1)*x(i)));
       dc1 = (exp(-l(1)^2*t(n))*A(1)*cos(l(1)*x(i))*l(1));
       d2c1 = -(exp(-l(1)^2*t(n))*A(1)*sin(l(1)*x(i))*l(1)^2);
       b(i,n+1)=b(i,n)+dt*(Db/Dc*(b(i+1,n)-2*b(i,n)+b(i-1,n))/dx^2+X0/Dc*(((b(i+1,n)-b(i-1,n))/(2*dx)/c1*dc1+b(i,n)/c1*d2c1-b(i,n)/c1^2*dc1^2)));
    end
    b(1,n+1)=0;
    b(N,n+1)=b(N-1,n+1);
end
hold on
for i=1:M
    plot(x,b(:,i))
end