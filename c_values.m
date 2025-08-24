close all
clear all
clc
L=1;
T=1;
N=101; %space
M=101; %time
% Db=linspace(10^-13,10^-12,100);
Db=10^-13;
Dc=10^-9;
X0=10^-11;
% X0=linspace(10^-11,10^-9,100);
dx=L/(N-1);
x=linspace(0,1,M);
t=linspace(0,1,N);
dt=T/(M-1);
l=[pi/2 3*pi/2 5*pi/2];
A=4.*(cos(l)-1)./(2.*l-sin(2.*l));
for n=1:M
    for i=1:N
       c1(n,i) = (1+exp(-l(1)^2*t(n))*A(1)*sin(l(1)*x(i)));
    end
end
k=0;
for i=[1 11 21 31 41 51 61 71 81 91 101]
k=k+1;
names{k}=sprintf('t = %ds',i);
plot(x,c1(i,:))
hold on
end
legend(names,'Location','best');
ylabel('Ammonia concentration')
xlabel('Length of the tube')
title('Distribution of ammonia along the length of the 1D pore')
figure(2)
surf(x,t,c1)
xlabel('x')
ylabel('t')
zlabel('c(x,t)')
grid on
title('Ammonia distribution');
view([150 25])
colorbar
colormap jet