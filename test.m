clear all
close all
clc
Dc=10e-5;
beta=0;
U=1/3*10^-4;
f0=1e-5;
M0=1;
gamma=U*(1-f0)/M0;
h=0.01;
zl=0;zh=100;
N=(zh-zl)/h+1;
z=zl:h:zh;
f_anly=sqrt(exp(M0.*z./6)./(2.*cosh(M0.*z./6)));
Rho_anly=(M0^2/6).*(exp(-M0.*z./12)./(2.*cosh(M0.*z./6)).^(1.5));
plot(z,f _anly)