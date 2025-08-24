N=101;
n=linspace(0,100,N);
dn=n(2)-n(1);
b0=linspace(0.2,0.4,N);
fun1 = @(b) root_sim_BACD(n,b,dn,N);
b = fsolve(fun1,b0);