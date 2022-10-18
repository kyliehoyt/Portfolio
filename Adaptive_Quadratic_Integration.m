clear all
close all
clc

a = -1;
b = 0.5;
fx = @(x)(x.^3 + x.^2 + 2);
[I,iter,seg] = Adaptive_Quad(fx,a,b)
integral(fx,a,b)

function [I,iter,seg] = Adaptive_Quad(f,a,b)
    iter = 0;
    seg = 2;
    err = 1.0;
    tol = 1e-8;
    intvecold = 0.0;
    % As you add more and more segments, the area will eventually stop
    % changing
    while mean(err)>tol
        iter = iter+1;
        step = (b-a)/seg;
        x = a:step:b;
        intvec = zeros(seg,1);
        for i=1:seg
            intvec(i) = 0.5*(x(i+1)-x(i))*(f(x(i))+f(x(i+1)));
        end
        intvec1 = reshape(intvec,[seg/2,2]);
        sum1 = sum(intvec1,2);
        err = (sum1-intvecold)/norm(sum1,2);
        intvecold = intvec;
        seg = seg*2;
    end
    I=sum(sum1);
    seg = seg/2;
end
