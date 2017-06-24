function [ x ] = denoise( y,tau )
c=size(y);
d1=rand(c);
d2=rand(c);
z=rand(c);
u=rand(c);
mu=1.5;
k=0;
%cost=zeros(300,1);
while(k<50)
x1=z+d1;
x2=u+d2;

x=1/2*(x1+x2);
z1=x-d1;
z=(mu*z1-1+((mu*z1-1).^2+4*mu*y).^0.5)/2/mu;
u1=x-d2;

u=tvd_mm(u1, tau/mu, 50);
%cost(k) = sum(x-y.*log(x))+ tau/2 * sum(abs(diff(x)));
d1=d1-x+z;
d2=d2-x+u;
k=k+1;
end
end

