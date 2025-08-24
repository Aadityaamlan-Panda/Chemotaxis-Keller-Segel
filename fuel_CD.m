function[c,f,s] = fuel_CD(r,t,u,dudr,beta,gamma)
c=[1;0;1];
f=[dudr(1)-u(1)*dudr(2);dudr(2);beta*dudr(3)];
s=[0;u(3)*u(1);-gamma*u(1)];
end