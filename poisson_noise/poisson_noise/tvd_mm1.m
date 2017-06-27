function [x] = tvd_mm1(y, lam, Nit)
% [x, cost] = tvd_mm(y, lam, Nit)
% Total variation denoising using majorization-minimization
% and banded linear systems.
%
% INPUT
%   y - noisy signal
%   lam - regularization parameter
%   Nit - number of iterations
%
% OUTPUT
%   x - denoised signal
%   cost - cost function history
%
% Reference
% 'On total-variation denoising: A new majorization-minimization
% algorithm and an experimental comparison with wavalet denoising.'
% M. Figueiredo, J. Bioucas-Dias, J. P. Oliveira, and R. D. Nowak.
% Proc. IEEE Int. Conf. Image Processing, 2006.

% Ivan Selesnick, selesi@nyu.edu, 2011
% Revised 2017
constant=1;
y = y(:);                                              % Make column vector
cost = zeros(1, Nit);                                  % Cost function history
N = length(y);

I = speye(N);
D = I(2:N, :) - I(1:N-1, :);
%DDT = D * D';

x = y;                                                 % Initialization
Dx = D*x;
%Dy = D*y;
c=1.1.*y.*log(x)+constant;
a=-y./x.^2+(c+y.*log(x))./x.^2;
i=find(abs(a)==-a);
    if(sum(i)~=0)
        a
        a(i)
        x(i)
        y(i)
        c(i)
        return;
    end
b=(y+x-2*(c+y.*log(x)))./x;
y_new=-b;
len=length(a);
a
%pause;
A=spdiags(a(:),0,len,len);

for k = 1:Nit
    F = sparse(1:N-1, 1:N-1, abs(Dx)/lam) + D*(A\D')/2;       % F : Sparse banded matrix
    x = A\y_new/2 - (A\D')/4*(F\D*(A\y_new));                 % Solve banded linear system
    %x
    %pause;
    c=1.1.*y.*log(x)+constant;
    a=-y./x.^2+(c+y.*log(x))./x.^2;
    i=find(abs(a)==-a);
    if(sum(i)~=0)
        a
        a(i)
        x(i)
        y(i)
        c(i)
        return;
    end
        
    b=(y+x-2*(c+y.*log(x)))./x;
    y_new=-b;
    %a
    %pause;
    A=spdiags(a(:),0,len,len);
    Dx = D*x;
    cost(k) = sum(x-y.*log(x)) + lam*sum(abs(Dx)); % cost function value
end
