clc;
clear all

X=60;
d=2*10^-6;
v=27*10^-6;
r=d/2;
t=linspace(0,50,51);
x=v*t;
m=0;

sol=pdepe(m,@lonepde,@loneic,@lonebc,x,t);

C=sol;

figure(3)
surf(x,t,C)
title('Attractant profile')
xlabel('x')
ylabel('t')
zlabel('C')
view([150 25])
colorbar;

th = 0:pi/50:2*pi;
disp=[0,0];
figure(4)
plot(disp(1),disp(2),'r*')
hold on
for time=1:49
    ratio=X/100*C(time,time)/C(time+1,time);
    num=rand;
    if ratio<num
        for i=v*time:5*r:v*(time+1)
        plot (r*cos(th)+disp(1),r*sin(th)+disp(2),'g-');
        ang=randsample(th,1);
        disp=disp+[r*cos(ang),r*sin(ang)];
        end
    else
        ang=randsample(th,1);
        plot ([disp(1),disp(1)+r*cos(ang)],[disp(2),disp(2)+r*sin(ang)],'b-')
        disp=disp+[r*cos(ang),r*sin(ang)];
    end
    plot(disp(1),disp(2),'ro')
end
title("Random Walk")
xlabel('x')
ylabel('y')
hold off



function [c,f,s] = lonepde(x,t,u,dudx)
D=8.0*10^-10;
c=1;
f=D*dudx;
s=0;
end

function u0 = loneic(x)
u0=0.001;
end

function [pl,ql,pr,qr] = lonebc(xl,ul,xr,ur,t)
C0=6938;
pl=ul-C0;
ql=0;
pr=ur;
qr=0;
end

