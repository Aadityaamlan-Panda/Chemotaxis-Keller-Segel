function [t1,b] = grid_variate_MoL(N,M,x,t,dx,Cch,T,d_ratio)
%% DEFINE CONSTANTS

% Db=linspace(10^-13,10^-12,100)
% X0=linspace(10^-11,10^-9,100)

Db= 2*10^-10;

Dc = 8.9*10^-10;

 % Manually enter d_ratio

k1 = 3.9*10^-10;
Ds=10^-9;

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

ci=0;

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

%% CALCULATING BACTERIAL CONCENTRATION

b0=ones(M, N).*1;
b0(1,:) = 1; 
b0(:,1) = 0;
tspan=[0 T];
   [t1,b]=ode45(@(t1,b1) bact_odefunc(t1,b1,N,d_ratio,dx,k1_by_Ds_by_Cch,k2_by_Cch,l,A,v,x), tspan, b0(1,:));
   %    b(1,:) = ones(size(b(1,:)));
% b(:,1) = zeros(size(b(:,1)));
end

