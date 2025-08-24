%% TIME IS ROW WHILE SPACE IS COLUMN IN M X N MATRIX DEFINED ALL THROUGHOUT
clc
clear all
close all

%% ENTER SPACE AND TIME DOMAIN

% MAXIMUM LENGTH
L=1;

% MAXIMUM TIME
T=2;

% NO OF POINTS IN SPACE DOMAIN
N=801;

% N is SPACE domain

% NO OF POINTS IN TIME DOMAIN
M=1001; 

% M is TIME domain

x=linspace(0,L,N);
t=linspace(0,T,M);

dt=T/(M-1);
dx=L/(N-1);

%% DEFINE CONSTANTS

% Db=linspace(10^-13,10^-12,100)
% X0=linspace(10^-11,10^-9,100)

Db= 2*10^-10;

Dc = 8.9*10^-10;



d_ratio = Db/Dc;

d_ratio = 0.001; % Manually enter d_ratio

k1 = 3.9*10^-10;
Ds=10^-9;

Cch=10^(-5);

k1_by_Ds_by_Cch= (k1/Ds)/Cch;

k2 = 5*10^-3;
k2_by_Cch = k2/Cch;

%% ENTER N NTH TERM/TERMS APPROXIMATION

v=100;

l=zeros(v,1);

for i=1:v
    l(i)=(2*i-1)*pi/2;
end

%% ENTER INITIAL CONCENTRATION IN ci

ci=0.01;

f=(1-ci)*4; %Initial factor

A=f.*(cos(l)-1)./(2.*l-sin(2.*l));

%% ENTER THE LENGTHS AT WHICH YOU WANT TO PLOT CONC VS TIME

l_points=[0.1 0.5 0.9];

%% ENTER THE TIMES AT WHICH YOU WANT TO PLOT BACTERIA VS LENGTH

t_points=[0.5 1 2];

%% CALCULATING C AND DC/DX AND 1/C*DC/DX

c=ones(M,N);
dc=zeros(M,N);
vc=zeros(M,N);
c1=ones(M,N);

for n=1:M
    for i=1:N
       c1(n,i) = c1(n,i) + exp(-l(1)^2*t(n))*A(1)*sin(l(1)*x(i));
       for k=1:v
       c(n,i) = c(n,i) + exp(-l(k)^2*t(n))*A(k)*sin(l(k)*x(i));
       dc(n,i) = dc(n,i) + exp(-l(k)^2*t(n))*A(k)*cos(l(k)*x(i))*l(k);
       end
       vc(n,i) = -k1_by_Ds_by_Cch*dc(n,i)/(c(n,i)+k2_by_Cch)^2;
    end
end

%% PLOTTING CONCENTRATION VS TIME AT DIFFERENT LENGTHS

k=0;
o=0;
hold on
for i=l_points./L.*(N-1)
    k=k+1;
    o=o+1;
    plot(t,c(:,i+1),"LineWidth",2)
    names{k}=sprintf('x = %0.2f*L_c',l_points(o));
    k=k+1;
    plot(t,c1(:,i+1),'--',"LineWidth",2)
    names{k}=sprintf('x = %0.2f*L_c (FTA)',l_points(o));
end
legend(names,'Location','best')
ylabel('Ammonia concentration')
xlabel('Time')
title('Distribution of ammonia with time')
hold off

%% CALCULATING BACTERIAL CONCENTRATION

b0=ones(M, N).*1;
b0(1,:) = 1; 
b0(:,1) = 0;
tspan=[0 T];
   [t1,b]=ode89(@(t1,b1) bact_odefunc(t1,b1,N,d_ratio,dx,k1_by_Ds_by_Cch,k2_by_Cch,l,A,v,x), tspan, b0(1,:));
   %    b(1,:) = ones(size(b(1,:)));
% b(:,1) = zeros(size(b(:,1)));
b(:,N)=b(:,N-1);
%% PLOTTING BACTERIAL CONCENTRATION

figure(2)
hold on
k=0;
for i=t_points./T.*(size(t1,1)-1)
    k=k+1;
    plot(x,b(floor(i)+1,:),"LineWidth",2)
    names1{k}=sprintf('t = %0.2f*t_d',t_points(k));
end
legend(names1,'Location','best')
ylabel('Bacteria concentration')
xlabel('Length of the tube')
title(sprintf('d ratio = %0.2s and Cch = 1.00e%0.1f',d_ratio, log10(Cch)))
hold off

%% PLOTTING VELOCITY

figure(3)
hold on
k=0;
for i=t_points./T.*(M-1)
    k=k+1;
    plot(x,vc(i+1,:),'LineWidth',2)
    names2{k}=sprintf('t = %0.2f*t_d',t_points(k));
end
legend(names2,'Location','best')
ylabel('Velocity')
xlabel('Length of the tube')
title('Velocity along the length of the 1D pore')
hold off

k=0;
figure(4)
hold on
for i=l_points./L.*(N-1)
    k=k+1;
    plot(t,vc(:,i+1),'LineWidth',2)
    names3{k}=sprintf('x = %0.2f*L_c',l_points(k));
end
legend(names3,'Location','best')
ylabel('Velocity')
xlabel('Time')
title('Velocity with time')
hold off

%% Quantitative Upwind Algorithm Check
c_=randn(M,N);
c_(1, :)=0; 
c_(:,1)=1;
for n = 1:M-1
    for i=2:N-1
        c_(n+1,i)=c_(n,i)+dt*((c_(n,i-1)-2*c_(n,i)+c_(n,i+1))/dx^2);
    end
    c_(n+1,N)=c_(n+1,N-1);
end

figure(5)
k=0;
o=0;
hold on
for i=l_points./L.*(N-1)
    k=k+1;
    o=o+1;
    plot(t,c(:,i+1),"LineWidth",2)
    names4{k}=sprintf('x = %0.2f*L_c Analytical',l_points(o));
    k=k+1;

    plot(t,c_(:,i+1),'--',"LineWidth",2)
    names4{k}=sprintf('x = %0.2f*L_c Algorithm',l_points(o));
end
legend(names4,'Location','best')
ylabel('Ammonia concentration')
xlabel('Time')
title('Distribution of ammonia with time')
hold off