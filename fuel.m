clear all
close all
clc
Dc=10e-5;
beta=0;
U=1/3*10^-4;
f0=0;
M0=1;
gamma=U*(1-f0)/M0;
m=0;
r=linspace(0,1,100);
t=linspace(0,1,100);
eqn= @(r,t,u,dudr)fuel_CD(r,t,u,dudr,beta,gamma);
ic=@fuel_CD_ic;
bc = @(rl,ul,rr,ur,t)fuel_CD_bc(rl,ul,rr,ur,t,f0,gamma,M0);
sol=pdepe(m,eqn,ic,bc,r,t);
bact_conc=sol(:,:,1);
chem_conc=sol(:,:,2);
fuel_conc=sol(:,:,3);
for i=[1 11 51 100]
    plot(r,bact_conc(i,:))
    hold on
end
legend('1','11','51','100');
hold off