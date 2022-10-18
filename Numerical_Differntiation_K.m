%%%% Assignment 2, Problem 1 by Kylie Hoyt %%%%

clear
clc
close all

c0 = 40.0;
ccrit = 1.0;

k = -0.05;
deltat = 1.0;

t(1) = 0.0;
c(1) = c0;
i=1;

while abs(c(i)-ccrit) >= 0.1
    i= i + 1;
    c(i) = c(i-1) + k*deltat*(c(i-1)^1.5+c(i-1));
    t(i) = t(i-1) + deltat;
end

plot(t,c,'ok');
hold on
xlabel('Time')
ylabel('Feed Mass')
disp('---------------------------Problem 1---------------------------')
fprintf('Using delta t = %1.1f, the feed mass in the reactor will be %3.2f at time t = %3.2f',deltat, c(i),t(i))

