function [pl,ql,pr,qr] = fuel_CD_bc(rl,ul,rr,ur,t,f0,gamma,M0)
pl=[ul(1);-M0/(1-f0)*(gamma+(1+f0+f0^2)/6-ul(3)^2/2);ul(3)-f0];
ql=[0;1;0];
pr=[ur(1);ur(2);ur(3)-1];
qr=[0;0;0];
end