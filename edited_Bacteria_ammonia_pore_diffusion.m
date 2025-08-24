clear all
clc


m=0;
x=linspace(0,1,1000);
t=linspace(0,1,1000);

sol=pdepe(m,@coupledpde,@coupledic,@coupledbc,x,t);

b=sol(:,:,1);

figure(1)
surf(xdom,tdom,b)
xlabel('x')
ylabel('t')
zlabel('b(x,t)')
grid on
title('Bacterial distribution');
view([150 25])
colorbar;

for j=[1 101 1000]
 ylim([0,1.2])
figure(2)
 hold on
 plot(x,b(j,:),LineWidth=2)
 title(["Particle Concentration Profile at t= ",num2str(t(j)),'*td'])
 xlabel("Characteristic length")
 ylabel('bacteria Concentration')
 
end
function [c,f,s]= coupledpde(x,t,u, dudx)
% Db=linspace(10^-13,10^-12,100);
Db=10^-12;
Dc=10^-9;
X0=10^-9;
% X0=linspace(10^-11,10^-9,100);
c=1;
s=0;
l=[pi/2 3*pi/2 5*pi/2];
b=2./l;
C0=0.9;
mob_ratio=X0/Dc;
d_ratio=Db/Dc;
f=(d_ratio)*dudx-u*(mob_ratio).*(l(1).*(1-C0).*b(1)*cos(l(1).*x).* exp(-1.*(l(1).^2).*t)./(((C0 +(1-C0).*(b(1)*sin(l(1).*x).* exp(-1.*(l(1).^2).*t))))));
end

function u0 = coupledic(x)
u0=1;
end
function [pl,ql,pr,qr]= coupledbc(xl,ul,xr,ur,t)
pl=ul;
ql=0;
pr=0;
qr=1;
end