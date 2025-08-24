clear all
close all
clc
%% ENTER SPACE AND TIME DOMAIN

% MAXIMUM LENGTH
L=1;

% MAXIMUM TIME
T=2;

% NO OF POINTS IN SPACE DOMAIN
N=101;
N1=(N-1)/2+1;

% N is SPACE domain

% NO OF POINTS IN TIME DOMAIN
M=1001; 
M1=2*M-1;

% M is TIME domain

x=linspace(0,L,N);
x1=linspace(0,L,N1);
t=linspace(0,T,M);
t1=linspace(0,T,M1);

dt=T/(M-1);
dt1=T/(M1-1);
dx=L/(N-1);
dx1=L/(N1-1);

t_points=[0.5 1 2];

Cch=10^-(5.8);
%% PLOTTING BACTERIAL CONCENTRATION
b=grid_variate(N,M,x,t,dt,dx,Cch);
b1=grid_variate(N,M1,x,t1,dt1,dx,Cch);
b2=grid_variate(N1,M,x1,t,dt,dx1,Cch);

figure(1)
hold on
k=0;
o=0;
for i=t_points./T.*(M-1)
    k=k+1;
    o=o+1;
    plot(x,b(i+1,:),"LineWidth",2)
    names1{k}=sprintf('t = %0.2f*t_d for original points',t_points(o));
    k=k+1;
    plot(x,b1(2*i+1,:),'--',"LineWidth",2)
    names1{k}=sprintf('t = %0.2f*t_d for double points in time',t_points(o));
end
legend(names1,'Location','best')
ylabel('Bacteria concentration')
xlabel('Length of the tube')
title('Grid independency check by doubling points in time')
hold off

figure(2)
hold on
k=0;
o=0;
for i=t_points./T.*(M-1)
    k=k+1;
    o=o+1;
    plot(x,b(i+1,:),"LineWidth",2)
    names1{k}=sprintf('t = %0.2f*t_d for original points',t_points(o));
    k=k+1;
    plot(x1,b2(i+1,:),'--',"LineWidth",2)
    names1{k}=sprintf('t = %0.2f*t_d for half points in space',t_points(o));
end
legend(names1,'Location','best')
ylabel('Bacteria concentration')
xlabel('Length of the tube')
title('Grid independency check for half points in space')
hold off

disp(min(min(b)))