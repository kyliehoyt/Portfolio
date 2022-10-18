% False Position method of finding roots by Kylie Hoyt
% A bracketing method

clear
clc
close all

fx = @(x)(cos(0.8*x) + 0.01*x^2);

figure(1)
x = -4:0.1:8;
plot(x,fx(x),'-r')
grid on

maxiter = 100;
tol = 0.01;

xl = -10;
delta = 1;
xr =  xl+delta;

% Find the bracket
while fx(xl)*fx(xr)>0
    xr = xr + delta;
end

% Find the root using False Position
err = 1;
iter = 0;
xmold = 0;
while err>tol && iter<maxiter
    iter = iter +1;
    xm = xr-(fx(xr)*(xr-xl))/(fx(xr)-fx(xl));
    xr = xm;
    err = abs((xm-xmold)/xm)*100;
    xmold = xm;
end

if (iter<=maxiter)
    fprintf('\n The root is %10.5f, iter = %d\n',xm,iter)
else
    fprint('\n Maximum iteration of %d exceeded, no root',maxiter)
end