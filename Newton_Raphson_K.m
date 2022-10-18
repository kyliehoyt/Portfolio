% Newton-Raphson method of finding roots by Kylie Hoyt
% An open method

clear all
clc
close all

fx = @(x)(2*x.^2+5*x-7);
df = @(xm,xmold)(fx(xmold)/(xmold-xm));

figure(1)
x = -3:0.1:5;
plot(x,fx(x),'-r');
grid on;
legend('f(x)');

tol = 0.01;
err = 1.0;
maxiter = 100;

iter = 0;
xmold = 5;
xm = -1;
dfm = 4;

while err>tol && iter<maxiter
    iter = iter + 1;
    xm = xmold-fx(xmold)/dfm;
    dfm = df(xm,xmold);
    err = abs((xm-xmold)/xm)*100;
    xmold = xm;
end
if abs(fx(xm)) <=0.05
    fprintf('Root = %3.5f',xm)
else
    disp('Root not found')
end
