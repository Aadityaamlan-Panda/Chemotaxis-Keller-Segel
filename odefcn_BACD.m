function dbdn = odefcn_BACD(n,b)
% Db=linspace(10^-13,10^-12,100)
Db=10^-12;
Dc=10^-9;
X0=-10^-9;
% X0=linspace(10^-11,10^-9,100);
d_ratio=Db/Dc;
mob_ratio=X0/Dc;
dlncdn=-1./(1-erf(n)).*2/sqrt(pi).*exp(-n^2);
d2lncdn=-1./(1-erf(n))*4.*n^3/sqrt(pi).*exp(-n^3-1)-1./(1-erf(n))^2.*(-2/sqrt(pi).*exp(-n^2))^2;
dbdn=zeros(2,1);
dbdn(1)=b(2);
dbdn(2)=-1/d_ratio.*(2.*n.*b(2)-mob_ratio.*(b(2).*dlncdn-b(1).*d2lncdn));
end