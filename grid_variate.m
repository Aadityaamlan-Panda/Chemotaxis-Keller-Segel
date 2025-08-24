function b = grid_variate(N,M,x,t,dt,dx,Cch)
%% DEFINE CONSTANTS

% Db=linspace(10^-13,10^-12,100)
% X0=linspace(10^-11,10^-9,100)

Db= 2*10^-10;

Dc = 8.9*10^-10;

d_ratio = Db/Dc;

d_ratio = 0.001; % Manually enter d_ratio

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

b=zeros(M, N);
b(1,:) = 1; 
b(:,1) = 0;
for n = 1:M-1
    for i=3:N-1
        b(n+1,i)=b(n,i)+dt*(-(b(n,i)*(3*vc(n,i)-4*vc(n,i-1)+vc(n,i-2))/(2*dx)+vc(n,i)*(3*b(n,i)-4*b(n,i-1)+b(n,i-2))/dx)+d_ratio*(b(n,i-1)-2*b(n,i)+b(n,i+1))/dx^2);
    end
    b(:,N)=b(:,N-1);
end