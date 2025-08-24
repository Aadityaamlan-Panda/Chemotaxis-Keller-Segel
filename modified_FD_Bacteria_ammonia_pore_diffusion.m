clear all
close all
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
b=zeros(N,M);
b(:,1)=1;
l=[pi/2 3*pi/2 5*pi/2];
A=4.*(cos(l)-1)./(2.*l-sin(2.*l));
if Db * dt / dx^2 > 0.5
    error('The time step is too large for stability.')
end
for n=1:M
    for i=1:N
       c1(i,n) = (1+exp(-l(1)^2*t(n))*A(1)*sin(l(1)*x(i)));
       dc1(i,n) = (exp(-l(1)^2*t(n))*A(1)*cos(l(1)*x(i))*l(1));
       vc(i,n) = dc1(i,n)/c1(i,n);
    end
end
for n=1:M-1
    for i=2:N-1
       b(i,n+1)=b(i,n)+dt*(Db/Dc*(b(i+1,n)-2*b(i,n)+b(i-1,n))/dx^2-X0/Dc*((b(i+1,n)*vc(i+1,n)-b(i,n)*vc(i,n))/(dx)));
    end
    b(1,n+1)=0;
    b(N,n+1)=b(N-1,n+1);
end
hold on
k=0;
for i=[1 11 101]
    k=k+1;
    plot(x,b(:,i))
    names{k}=sprintf('t = %ds',i);
    xlim([0 1]);
    ylim([0 2])
end
legend(names,'Location','best')
ylabel('Bacteria concentration')
xlabel('Length of the tube')
title('Distribution of bacteria along the length of the 1D pore')
figure(2)
surf(x,t,b)
xlabel('x')
ylabel('t')
zlabel('b(x,t)')
grid on
title('Bacterial distribution');
view([150 25])
colorbar
colormap jet