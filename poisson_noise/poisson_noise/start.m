
x=load('blocks.txt');
x=x+2*abs(min(x));
%testing signal below and above
%a=[1 1 1 1 1 1 1 1 1 1  1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 4 4 4 4 4 4 4 4 4 4 4 4 4 4  4 4 4 4 4 4 4 4 4 4 ];
%x=(repmat(a,1,4))';
epsilon=0.00001;
N=256;
Nit=200;
figure(1);
subplot(2,1,1);
plot(x);
title('noise-free signal');
xlabel('time');
%poisson distribution for noise
y=noise(x)+epsilon;
subplot(2,1,2);
plot(y);
title('noisy signal');
xlabel('time');

    
lam = 1.6;                         % lam: regularization parameter
%denoising
mu=lam/50;
[z]=denoise(y,lam,Nit,mu);
figure;
subplot(2,1,1);
plot(z);
title('denoised(using admm) signal');
xlabel('time');
%comparison with alogo for  gaussian noise
z = tvd_mm(y,lam,Nit);
subplot(2,1,2);
plot(z);
title('denoised(using MM for gaussian) signal');
ylabel('time');


