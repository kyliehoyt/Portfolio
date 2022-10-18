clear all
close all
clc

a = -1;
b = 0.5;
fx = @(x)(x.^3 + x.^2 + 2);
Composite_Trap(fx,a,b)
integral(fx,a,b)
xx = a:0.05:b;
plot(xx,fx(xx),'-k')

function [sum] = Composite_Trap(f,a,b)
    seg = 2;
    err = 1.0;
    tol = 1e-8;
    sumold = 0.0;
    % As you add more and more segments, the area will eventually stop
    % changing
    while err>tol
        step = (b-a)/seg;
        x = a:step:b;
        sum = 0;
        for i=1:seg
            sum = sum + (x(i+1)-x(i))*(f(x(i))+f(x(i+1)))/2;
        end
        err = abs((sum-sumold)/sum);
        sumold = sum;
        seg = seg*2;
    end
end