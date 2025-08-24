function [c,f,s]= coupledpde_BACD(x,t,u, dudx)
Db=10^-12;
Dc=10^-9;
X0=10^-9;
% X0=linspace(10^-11,10^-9,100);
% Db=linspace(10^-13,10^-12,100);
c=1;
s=0;
l=[pi/2 3*pi/2 5*pi/2];
A=4.*(cos(l)-1)./(2.*l-sin(2.*l));
mob_ratio=-1;
d_ratio=0.01;
f=(d_ratio).*dudx-u*(mob_ratio)./(1+exp(-l(1)^2.*t).*A(1).*sin(l(1).*x)).*(exp(-l(1)^2.*t).*A(1).*cos(l(1).*x).*l(1));
end