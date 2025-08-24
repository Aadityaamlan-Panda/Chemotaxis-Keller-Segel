clear all
clc
close all

m=0;
x=linspace(0,1,1001);
t=linspace(0,1,1001);

sol=pdepe(m,@coupledpde_BACD,@coupledic_BACD,@coupledbc_BACD,x,t,options);

b=sol(:,:,1);

figure(1)
surf(x,t,b)
xlabel('x')
ylabel('t')
zlabel('b(x,t)')
grid on
title('Bacterial distribution');
view([150 25])
colorbar
colormap jet

for j=[1 11 101]
 ylim([0,1.2])
 figure(2)
 hold on
 plot(x,b(j,:),LineWidth=2)
 title("Particle Concentration Profile at t= (1, 11, 101)*td")
 xlabel("Characteristic length")
 ylabel('bacteria Concentration')
end