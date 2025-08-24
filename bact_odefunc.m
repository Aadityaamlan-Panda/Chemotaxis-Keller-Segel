function dbdt = bact_odefunc(t,b,N,d_ratio,dx,k1_by_Ds_by_Cch,k2_by_Cch,l,A,v,x)
dbdt = zeros(N, 1);
c=ones(N,1);
dc=zeros(N,1);
for i = 1:N
       for k=1:v
       c(i) = c(i) + exp(-l(k)^2.*t).*A(k).*sin(l(k)*x(i));
       dc(i) = dc(i) + exp(-l(k)^2.*t).*A(k).*cos(l(k)*x(i)).*l(k);
       end
       vc1(i) = -k1_by_Ds_by_Cch*dc(i)/(c(i)+k2_by_Cch)^2;
end
for i = 2:N-1
    if i==2
        dbdt(i)=-(b(i)*(vc1(i+1)-vc1(i-1))/(2*dx)+vc1(i)*(-b(i-1)+b(i+1))/(2*dx))+d_ratio*(b(i-1)-2*b(i)+b(i+1))/dx^2;
    else
         dbdt(i)=-(b(i)*(3*vc1(i)-4*vc1(i-1)+vc1(i-2))/(2*dx)+vc1(i)*(3*b(i)-4*b(i-1)+b(i-2))/(2*dx))+d_ratio*(b(i-1)-2*b(i)+b(i+1))/dx^2;
    end
end
disp("done")
end
