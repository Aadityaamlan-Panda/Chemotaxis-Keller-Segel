function eq=myfunc(c_,M,N,dx,dt)
eq=zeros(M,N);
eq(1,:) = c_(1, :); 
eq(:,1) = c_(:,1)-1;
for n = 1:M-1
    for i=2:N-1
        %b(n+1,i)=b(n,i)+dt*(-mob_ratio*(b(n,i)*(vc(n,i+1)-vc(n,i))/dx+vc(n,i)*(b(n,i+1)-b(n,i))/dx)+d_ratio*(b(n,i-1)-2*b(n,i)+b(n,i+1))/dx^2);
        eq(n+1,i)=-c_(n+1,i)+c_(n,i)+dt*((c_(n,i-1)-2*c_(n,i)+c_(n,i+1))/dx^2);
    end
    eq(n+1,N)=c_(n+1,N)-c_(n+1,N-1);
end