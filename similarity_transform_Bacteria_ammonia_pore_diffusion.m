clc
clear all
close all
Dc = 10^-9;
b0 = [0;1e-5];
nspan = [0 100];
[n,b] = ode45(@odefcn_BACD, nspan, b0);
c=1-erf(n);
hold on
k=0;
for i=1:100:501
    plot(n,b(:,1))
    k=k+1;
    names{k}=sprintf('t = %ds',i);
end
legend(names,'Location','best')
ylabel('Bacteria concentration')
xlabel('eta')
title('Distribution of bacteria with eta')
% tq = linspace(0.1, 20*10^9, 501);
% k=0;
% for i=1:length(tq)
%     for j=1:length(n)
%         u(i,j)=n(j)*sqrt(4*Dc*tq(i));
%     end
% end
% hold on
% grid on
% b(1,1)=1;
% for i=1:100:501
%     plot (u(i,:),b(:,1))
%     k=k+1;
%     names{k}=sprintf('t = %ds',i);
%     xlim([0 55]);
% end
% legend(names,'Location','best')
% ylabel('Bacteria concentration')
% xlabel('Length of the tube')
% title('Distribution of bacteria along the length of the 1D pore')
% hold off
% figure(2)
% hold on
% grid on
% k=0;
% for i=1:100:501
%     plot(u(i,:),c)
%     k=k+1;
%     names{k}=sprintf('t = %ds',i);
%     xlim([0 55]);
% end
% legend(names,'Location','best')
% ylabel('Solute Concentration')
% xlabel('Length of the tube')
% title('Diffusion of ammonia along the length of the 1D pore')