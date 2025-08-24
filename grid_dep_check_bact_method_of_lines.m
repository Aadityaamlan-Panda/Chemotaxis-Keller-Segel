clear all
close all
clc
%% ENTER SPACE AND TIME DOMAIN

% MAXIMUM LENGTH
L=1;

% MAXIMUM TIME
T=2;

% NO OF POINTS IN SPACE DOMAIN
N=401;
N1=801;

% N is SPACE domain

% NO OF POINTS IN TIME DOMAIN
M=1001; 


% M is TIME domain

x=linspace(0,L,N);
x1=linspace(0,L,N1);
t=linspace(0,T,M);

dt=T/(M-1);

dx=L/(N-1);
dx1=L/(N1-1);

t_points=[0.5 1 2];

Cch=10^(-5);
d_ratio=0.001;
%% PLOTTING BACTERIAL CONCENTRATION
[t1,b]=grid_variate_MoL(N,M,x,t,dx,Cch,T,d_ratio);
[t2,b1]=grid_variate_MoL(N1,M,x1,t,dx1,Cch,T,d_ratio);

figure(1)
hold on
k=0;
for i=t_points./T.*(size(t1,1)-1)
    k=k+1;
    plot(x,b(floor(i)+1,:),"LineWidth",2)
    names1{k}=sprintf('t = %0.2f*t_d for %d points',t_points(k),N);
end

for i=t_points./T.*(size(t2,1)-1)
    plot(x1,b1(floor(i)+1,:),'--',"LineWidth",2)
    k=k+1;
    names1{k}=sprintf('t = %0.2f*t_d for %d points',t_points(k-size(t_points,2)),N1);
end
legend(names1,'Location','best')
ylabel('Bacteria concentration')
xlabel('Length of the tube')
title(sprintf('d ratio = %0.2s and Cch = 1.00e%0.1f',d_ratio, log10(Cch)))
hold off