clc;
clear all
t=linspace(0,10,101);
x=linspace(0,1,101);
m=0;

sol=pdepe(m,@coupledpde,@coupledic,@coupledbc,x,t);

u_=sol(:,:,1);
v_=sol(:,:,2);

figure(1)
surf(x,t,u_)
xlabel('x')
ylabel('t')
zlabel('u_')
title("Bacterial density")
view([150 25])
colorbar;

figure(2)
surf(x,t,v_)
xlabel('x')
ylabel('t')
zlabel('v_')
title('Sensing molecule concentration')
view([150 25])
colorbar;

function [c,f,s] = coupledpde(x,t,u,dudx)
D1=0.1; D2=1; a=10; b=1;
ctc=0.002;
u1=1;
c=ones(2,1);
f=[D1*dudx(1)-ctc*u1*dudx(2);D2*dudx(2)];
s=[0;a*u(1)-b*u(2)];
end

function u0 = coupledic(x)
u0=[1;1+0.1.*exp(-10.*x.^2)];
end

function [pl,ql,pr,qr] = coupledbc(xl,ul,xr,ur,t)
pl=zeros(2,1);
ql=ones(2,1);
pr=zeros(2,1);
qr=ones(2,1);
end

