% Fixed Point method of finding roots by Kylie Hoyt
% An open method

clear all
clc
close all

fx = @(x)(cos(x));
gx = @(x)(fx(x) + x);

figure(1)
x = -3:0.1:5;
plot(x,fx(x),'-r')
hold on
plot(x,gx(x),'-g')
plot(x,x,'-k')
grid on;
legend('f(x)','g(x)','y = x')

tol = 0.01;
err = 1.0;
maxiter = 100;

iter = 0;
xmold = 0;
plot(xmold,gx(xmold),'bo')
while err>tol && iter<maxiter
    iter = iter + 1;
    xm = gx(xmold);
    xms(iter) = xm;
    err = abs((xm-xmold)/xm)*100;
    xmold = xm;
end
if abs(fx(xm)) <=1.0
    fprintf('Root = %3.5f',xm)
else
    disp('Root not found')
end
plot(xms,gx(xms),'bo')
plot(xms,xms,'bo')
