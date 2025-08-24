clear all
close all
clc
Dc=10e-5;
beta=0;
U=1/3*10^-4;
f0=0;
M0=1;
gamma=U*(1-f0)/M0;
h=0.001;
zl=0;zh=1000;
N=(zh-zl)/h+1;
z=zl:h:zh;
dRhodz=@(Rho,s,U) Rho*(s-U);
dsdz=@(Rho,f) -f*Rho;
dfdz=@(Rho,U,gamma) gamma*Rho/U;
s(1)=M0/(1-f0)*(gamma+(1+f0+f0^2)/6-f0^2/2);
Rho(1)=1e-5;
f(1)=f0;
for n=1:N-1
    Rho_k1=h*dRhodz(Rho(n),s(n),U);
    f_k1=h*dfdz(Rho(n),U, gamma);
    s_k1=h*dsdz(Rho(n),f(n));
    Rho_k2=h*dRhodz(Rho(n)+Rho_k1/2,s(n)+s_k1/2,U);
    f_k2=h*dfdz(Rho(n)+Rho_k1/2,U, gamma);
    s_k2=h*dsdz(Rho(n)+Rho_k1/2,f(n)+f_k1/2);
    Rho(n+1)=Rho(n)+Rho_k2;
    f(n+1)=f(n)+f_k2;
    s(n+1)=s(n)+s_k2;
   
end
Rho_anly=(U^2/(6*gamma^2)).*(f-f0).*(1-f).*(1+f+f0);
figure(1)
plot(z,Rho,'LineWidth',3)
hold on
plot(z,Rho_anly)
legend('\rho','\rho anly')
xlabel('\xi')
ylabel('\rho')
title(sprintf('Plot for bacterial density at zdiff = %0.2f and zmax = %d',h,zh))
hold off
figure(2)
plot(z,f,'LineWidth',3)
xlabel('\xi')
ylabel('f')
title('Plot for fuel')
figure(3)
plot(z,abs(Rho-Rho_anly)./Rho)
xlabel('\xi')
ylabel('\delta\rho')
title('Error')