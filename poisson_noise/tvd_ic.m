function [x,J] = tvd_ic(y,lam,Nit)
% [x,J] = tvd_ic(y,lam,a,Nit)
% Total variation denoising using
% iterative clipping algorithm.
% INPUT
%   y - noisy signal (row vector)
%   lam - regularization parameter
%   Nit - number of iterations
% OUTPUT
%   x - result of denoising
%   J - objective function

y = y(:)';              % row vector
J = zeros(1,Nit);       % objective function
N = length(y);
z = zeros(1,N-1);       % initialize z
% alpha = 4;
alpha = 3;
T = lam/2;
for k = 1:Nit
    x = y - [-z(1) -diff(z) z(end)];      % y - D' z
    J(k) = sum(abs(x-y).^2) + lam * sum(abs(diff(x)));
    z = z + 1/alpha * diff(x);            % z + 1/alpha D z
    z = max(min(z,T),-T);                 % clip(z,T)
end
