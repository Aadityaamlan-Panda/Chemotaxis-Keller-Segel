function eq1 = root_sim_BACD(n,b,dn,N)
% Db=linspace(10^-13,10^-12,100)
Db=10^-12;
Dc=10^-9;
X0=10^-9;
% X0=linspace(10^-11,10^-9,100);
d_ratio=Db/Dc;
mob_ratio=X0/Dc;
dlncdn=-1./(1-erf(n)).*2/sqrt(pi).*exp(-n.^2);
d2lncdn=-1./(1-erf(n))*4.*n.^3/sqrt(pi).*exp(-n.^3-1)-1./(1-erf(n)).^2.*(-2/sqrt(pi).*exp(-n.^2)).^2;
eq1=zeros(1,N);
eq1(1)=b(1);
eq1(N)=b(N);
for i=2:N-1
eq1(i)=(b(i-1)+b(i+1)-b(i))/dn^2+1/d_ratio*(2*n(i)*(b(i+1)-b(i-1))/(2*dn)-mob_ratio*((b(i+1)-b(i-1))/(2*dn)*dlncdn(i)-b(i)*d2lncdn(i)));
end
end