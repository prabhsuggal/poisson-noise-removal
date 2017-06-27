function [ x ] = denoise( y,tau,Nit,mu )
c=size(y);
d1=zeros(c);
d2=zeros(c);
z=y;
u=y;
%mu=tau/50;

cost=zeros(50,1);
for i=1:Nit
x1=z+d1;
x2=u+d2;

x=1/2*(x1+x2);                    % subproblem 1 solution due to admm 
z1=x-d1;
z=(mu*z1-1+((mu*z1-1).^2+4*mu*y).^0.5)/2/mu; % subproblem 2 solution due to admm 
u1=x-d2;
%u=(tvd_ic(u1,tau/mu,50))';          this was done for checking the better
                                        %b\w MM and iterative clipping
%size(u)
u=tvd_mm(u1, tau/mu, 50);                    % subproblem 3 solution due to admm 
cost(i)=sum(x-y.*log(abs(x)))+ tau * sum(abs(diff(x)));    %cost function
d1=d1-x+z;
%size(x)
%size(d2)
%size(u)
d2=d2-x+u;
end
figure;
plot(cost);
title('cost function');
ylabel('no. of iterations');
end

