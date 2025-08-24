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
dRhodz=@(Rho,s,U) Rho*(s-U);
dsdz=@(Rho,f) -f*Rho;
dfdz=@(Rho,U,gamma) gamma*Rho/U;
s(1)=M0/(1-f0)*(gamma+(1+f0+f0^2)/6-f0^2/2);
Rho(1)=1e-5;
f(1)=f0;
for n=1:N-1
    Rho_k1=dRhodz(Rho(n),s(n),U);
    Rho_k2=dRhodz(Rho(n)+h*Rho_k1/2,s(n),U);
    Rho_k3=dRhodz(Rho(n)+h*Rho_k2/2,s(n),U);
    Rho_k4=dRhodz(Rho(n)+h*Rho_k3,s(n),U);
    Rho(n+1)=Rho(n)+h/6*(Rho_k1+2*Rho_k2+2*Rho_k3+Rho_k4);

    s_k1=dsdz(Rho(n),f(n));
    s_k2=dsdz(Rho(n)+h*Rho_k1/2,f(n));
    s_k3=dsdz(Rho(n)+h*Rho_k2/2,f(n));
    s_k4=dsdz(Rho(n)+h*Rho_k3/2,f(n));
    s(n+1)=s(n)+h/6*(s_k1+2*s_k2+2*s_k3+s_k4);

    f_k1=dfdz(Rho(n),U, gamma);
    f_k2=dfdz(Rho(n)+h*Rho_k1/2,U, gamma);
    f_k3=dfdz(Rho(n)+h*Rho_k1/2,U, gamma);
    f_k4=dfdz(Rho(n)+h*Rho_k1/2,U, gamma);
    f(n+1)=f(n)+h/6*(f_k1+2*f_k2+2*f_k3+f_k4);
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