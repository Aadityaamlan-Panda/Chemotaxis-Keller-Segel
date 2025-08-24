function eq = root_BACD(d_ratio,mob_ratio,b,M,N,vc,dt,dx)
eq=zeros(M,N);
eq(1,:)=b(1,:)-1;
eq(:,1)=b(:,1);
for n=1:M-1
    for i=2:N-1
       eq(n+1,i)=-(b(n+1,i)-b(n,i))/dt+d_ratio*(b(n,i+1)-2*b(n,i)+b(n,i-1))/dx^2-mob_ratio*((b(n,i+1)*vc(n,i+1)-b(n,i)*vc(n,i))/(dx));
    end
    eq(n+1,N)=b(n+1,N)-b(n+1,N-1);
end
end