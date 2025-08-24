function eq = root3_BACD(d_ratio,b,M,N,vc,dt,dx)
eq=zeros(M,N);
eq(1,:)=b(1,:)-1;
eq(:,1)=b(:,1);
for n=1:M-1
for i = 2:N-1
    if i==2
        eq(n+1,i)=-(b(n+1,i)-b(n,i))/dt-(b(i)*(vc(i+1)-vc(i))/(dx)+vc(i)*(-b(i)+b(i+1))/dx)+d_ratio*(b(i-1)-2*b(i)+b(i+1))/dx^2;
    else
         eq(n+1,i)=-(b(n+1,i)-b(n,i))/dt-(b(i)*(3*vc(i)-4*vc(i-1)+vc(i-2))/(2*dx)+vc(i)*(3*b(i)-4*b(i-1)+b(i-2))/2*dx)+d_ratio*(b(i-1)-2*b(i)+b(i+1))/dx^2;
    end
end
    eq(n+1,N)=b(n+1,N)-b(n+1,N-1);
end
end
