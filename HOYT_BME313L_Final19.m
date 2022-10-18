%%%============== Final Exam =============%%%
% Kylie Hoyt

%% ------------Notes--------------- %%

% Total points: 100 (Assigned on 12/09/19: 9:00 AM)

% (1)STRICT DEADLINE: 9:00 AM Saturday (12/14/19).
% The deadline will be strictly enforced with a zero for late submission.

% (2) Questions
% Question will be entertained till 6:00 PM Wednesday (12/11/19).
% I will be holding office hours on Monday (12/09/19) and Tuesday (12/10/19)
% from 2 to 5 PM. Please use these to ask any question in person. Emails
% after 6:00 PM will be answered the next day.

% (3) Code neatness and formatted output
% The code must be written nicely with the function definitions at the
% bottom of your *.m script. The expected answere must be ouput on the
% screen where necessary.

%%
clear
clc
close all
fig = 0;
fprintf('BME313L Fall Final by Kylie Hoyt\n\n')
%--------------------------------------------------------------------------

%% Problem 1 (20 points +  5 for code neatness)

fprintf("---------------------------Problem 1a---------------------------\n")
% 1a. You are given the following data
% What polynomial order function will fit this data well? Use fprintf to 
% display your answer in the command window.
% Output the regressed polynomial coefficients, and r^2 value in the
% command window.
% Display your regressed curve in the plot as a red line with the original
% data as blue circles.

% GIVEN DATA
load P1a_data.mat
y1 = y1+5;

% POLYNOMIAL REGRESSION
fprintf('Data x1 and y1: 4th order polynomial:\n')
a1 = polyfit(x1,y1,4); % returns the coefficients to vector a
yapprox = a1(1)*x1.^4+a1(2)*x1.^3+a1(3)*x1.^2+a1(4)*x1+a1(5); % finds the yapprox using the coefficents from polyfit
fprintf('\tYapprox = %2.2fx^4 + %2.2fx^3 + %2.2fx^2 + %2.2fx + %2.2f\n',a1(1),a1(2),a1(3),a1(4),a1(5))

% R^2
Sr = sum((y1-yapprox).^2);
St = sum((y1-mean(y1)).^2);
rsq = 1-Sr/St;
fprintf('\tR^2 = %1.4f\n',rsq)

% PLOT
fig = fig+1;
figure(fig)
plot(x1,y1,'ob')
hold on
plot(x1,yapprox,'-r')
xlabel('x')
ylabel('y')
title('Problem 1a')

fprintf("\n---------------------------Problem 1b---------------------------\n")
% 1b. You are given the following multiple linear data 
% Output the regressed multiple linear coefficients, and r^2 value in the
% output window.

% GIVEN DATA
load P1b_data.mat

% MULTIPLE LINEAR REGRESSION
[a2,rsq2] = multiple_reg(y2,x2);
fprintf('Data x2 and y2: Multiple Linear Regression:\n')
fprintf('\tYapprox = %2.3f + %2.3f*x1 + %2.3f*x2\n',a2(1),a2(2),a2(3))
fprintf('\tR^2 = %1.4f\n',rsq2)

% PLOT
fig = fig+1;
figure(fig)
scatter3(x2(:,1),x2(:,2),y2)
xlabel('x1(1)')
ylabel('x1(2)')
zlabel('y1')
title('Problem 1b')

%--------------------------------------------------------------------------

%% Problem 2 (20 points + 5 for code neatness)

% GIVEN DATA
f2  = @(x)(1.*(x<=0.5)+0.*(x>0.5));
% bounded by [a=0,b=1] given by
a = 0;
b = 1;

fprintf("\n---------------------------Problem 2a---------------------------\n")
% 2a. What is the area under the curve using trapezoidal rule with 1,2,3,
% and 4  segements. Report the answer in the command window. 
% You will notice that the area alternates between the correct value 0.5 
% and different from 0.5. Explain why is this happening. Output the
% comment as an fprintf statement in the command window.

% TRUTH
fprintf('True integral = %1.4f\n',integral(f2,a,b));

% TRAPEZOIDAL AREAS
fprintf('Trapezoidal sum using %d segments: %1.4f\n',1,Segmented_Trap(f2,a,b,1))
fprintf('Trapezoidal sum using %d segments: %1.4f\n',2,Segmented_Trap(f2,a,b,2))
fprintf('Trapezoidal sum using %d segments: %1.4f\n',3,Segmented_Trap(f2,a,b,3))
fprintf('Trapezoidal sum using %d segments: %1.4f\n',4,Segmented_Trap(f2,a,b,4))

% ALTERNATING AREA
fprintf('\n\tThe area alternates between the correct value and an incorrect value\n')
fprintf('because of the spike and the location of the segments. When the spike\n')
fprintf('is in the middle of a segment, i.e. the number of segments is odd, the\n')
fprintf('error in the area evens out on either side of the spike. When there is\n')
fprintf('a data point at 0.5, i.e. the number of segments is even, the area to the\n')
fprintf('right of the spike is overestimated until the step size gets down to 0.01\n')
fprintf('which would be at 100 segments.\n')

fprintf("\n---------------------------Problem 2b---------------------------\n")
% 2b. Use the Trapezoidal and Simpson's rule to approximate the area under
% the curve with a error tolerance of 1.E-2. Report the number of
% segments to achieve this in the command window.

% COMPOSITE TRAPEZOIDAL
[Tsum,Tseg] = Composite_Trap(f2,a,b);
fprintf('\nTrapezoidal Rule:\n')
fprintf('\tSum: %1.4f\n \tSegments: %d\n',Tsum,Tseg)

% COMPOSITE SIMPSONS
[Ssum,Sseg] = Composite_Simpsons(f2,a,b);
fprintf('\nSimpsons Rule:\n')
fprintf('\tSum: %1.4f\n \tSegments: %d\n',Ssum,Sseg)

% PLOT
x3 = a:0.01:b;
y3 = f2(x3);
fig = fig+1;
figure(fig)
plot(x3,y3,'-b')
xlabel('x2')
ylabel('y2')
title('Problem 2')
ylim([-0.1,2])

%--------------------------------------------------------------------------

%% Problem 3 (30 points + 5 for code neatness)

% GIVEN DATA
load P3_data.mat

fprintf("\n---------------------------Problem 3a---------------------------\n")
% 3a. Calculate the forward, backward, and centered difference
% approximations for the first derivative df/dx (or dy/dx) at  appropriate
% points. Output the result in the command window as a table. 

% DISCRETE DIFFERENCES
fprintf('Discrete Points: h = 1\n')
f_df = forward_df(x3,y3);
b_df = backward_df(x3,y3);
c_df = centered_df(x3,y3);

% TABLE
Data = [f_df; b_df; c_df];
x1df = Data(:,1);
x2df = Data(:,2);
x3df = Data(:,3);
x4df = Data(:,4);
Cols = ["x_0";"x_1";"x_2";"x_3"];
Rows = ["Forward","Backward","Centered"];
Discrete_DF = table(x1df,x2df,x3df,x4df,'VariableNames',Cols,'RowNames',Rows)

fprintf("---------------------------Problem 3b---------------------------\n")
% 3b. Since (x3,y3) dataset consists of four points, we can fit a cubic
% polynomial to this dataset. Use Lagrange interpolation to get a 
% third order/cubic interpolant so that you can calculate a better 
% approximation of the derivative using a smaller increments 'h' on the 
% interpolated cubic function. 
% Report the interpolant coefficients in the command window. Since you now
% have the fitting function evaluate the foward, backward, and central
% difference with an 'h' of 0.01 at x = 2 and x = 2.5 Report the answer in the command window.

% LAGRANGE FUNCTION AND COEFFICIENTS
fprintf('Data x3 and y3: Lagrange Interpolation:\n')
[Lf,Lc] = Lagrange(x3,y3);
fprintf('\tLagrangeFunc = %2.1fx^3 + %2.1fx^2 + %2.1fx + %2.1f\n',Lc(1),Lc(2),Lc(3),Lc(4))

% LAGRANGE CONTINUOUS DIFFERENTIATION
fprintf('\nContinuous Function: h = 0.01\n')
[Lagrange_DF] = Continuous_Diff([2, 2.5],Lf,0.01)

% PLOT
fig = fig+1;
figure(fig)
plot(x3,y3,'ob')
hold on
x = min(x3):0.01:max(x3);
plot(x,Lf(x),'-r')
xlabel('x3')
ylabel('y3')
title('Problem 3')

%--------------------------------------------------------------------------

%% Problem 4 (10 points + 5 for code neatness)

fprintf("---------------------------Problem 4a---------------------------\n")
% 4a. Which type of linear solver (direct or iterative) will work here and Why?
% Use the appropriate solver and output the answer in the command window.

% GIVEN DATA
mat_A = [5 1 0; 2 9 3; 4 2 1];
b = [1 1 1]';

% GAUSSIAN ELIMINATION
fprintf('Direct Solver: Gaussian Elimination:\n')
[x] = Gaussian_Elim(mat_A,b)
Truth = mat_A\b
fprintf('\tOnly a direct solver would work here because matrix A is not diagonally\n')
fprintf('dominant which is a prerequisite for iterative solvers like the Jacobi Method.\n')

fprintf("\n---------------------------Problem 4b---------------------------\n")
% 4b. What would be the values of the regressed coefficents? Why?
% Hint: You do not need to solve anything, you only need to know the link
% between multiple linear regression and polynomial interpolation.
% Remember the construction of X in X'*X*(coefficient_vector) = X'*y in the
% notes.

% GIVEN DATA
a = 1;
b = 1;
c = 2;
f4 = @(x)(a*x.^2+b*x+c);
% We generate data point pairs (x4,y4) using,
x4 = (-1:0.1:1)';
y4 = f4(x4);
% If we used multiple linear regression using,
x4_new = [x4 x4.^2];
y4_new = y4;

% ANSWER
fprintf('\tThe coefficients would be 1, 1, and 2 exactly because polynomial\n')
fprintf('interpolation and multiple linear regression work in the same way.\n')
fprintf('Multiple linear regression works under the assumption that columns\n')
fprintf('are not linearly related. In this case they are related by the power\n')
fprintf('of 2 so multiple linear regression can still be used for the polynomial\n')
fprintf('function and it will regress to the actual coefficients.\n')


%% FUNCTIONS
%% Problem 1 Functions
%---------------------------Multiple Linear Regression
function [a,rsq] = multiple_reg(y,x)
    Z = [ones(length(x),1) x];
    a = (Z'*Z)\(Z'*y);
    Sr = sum((y-Z*a).^2);
    St = sum((y-mean(y)).^2);
    rsq = 1-Sr/St;
end

%% Problem 2 Functions
%---------------------------Trapezoidal Sum w/ Select Segmentation
function [sum] = Segmented_Trap(f,a,b,seg)
    step = (b-a)/seg;
    x = a:step:b;
    sum = 0;
    for i=1:seg
        sum = sum + 0.5*(x(i+1)-x(i))*(f(x(i))+f(x(i+1)));
    end
end

%---------------------------Composite Trapezoidal Sum
function [sum,seg] = Composite_Trap(f,a,b)
    seg = 0;
    err = 100;
    tol = 1e-2;
    sumold = 0.0;
    % As you add more and more segments, the area will eventually stop
    % changing
    while err>tol
        seg = seg + 1;
        step = (b-a)/seg;
        x = a:step:b;
        sum = 0;
        for i=1:seg
            sum = sum + 0.5*(x(i+1)-x(i))*(f(x(i))+f(x(i+1)));
        end
        err = abs((sum-sumold)/sum);
        sumold = sum;
    end
end

%---------------------------Composite Simpson's Rule
function [sum,seg] = Composite_Simpsons(f,a,b)
    seg = 1;
    err = 100;
    tol = 1e-2;
    sumold = 0.0; 
    while err>tol
        seg = seg + 1;
        step = (b-a)/seg;
        x = a:step:b;
        weights = step*[1/6 2/3 1/6];
        sum = 0;
        for i=1:seg
            x_seg = [x(i) 0.5*(x(i)+x(i+1)) x(i+1)];        
            sum = sum + weights*f(x_seg)';  
        end
        err = abs((sum-sumold)/sum);
        sumold = sum;
    end
end

%% Problem 3 Functions
%---------------------------Forward Difference
function [df] = forward_df(x,y)
    df = zeros(1,length(x));
    for i = 1:length(x)-1
        h = x(i+1)-x(i);
        df(i) = (y(i+1)-y(i))/h;
    end
end

%---------------------------Backward Difference
function [df] = backward_df(x,y)
    df = zeros(1,length(x));
    for i = 2:length(x)
        h = x(i-1)-x(i);
        df(i) = (y(i-1)-y(i))/h;
    end
end

%---------------------------Centered Difference
function [df] = centered_df(x,y)
    df = zeros(1,length(x));
    for i = 2:length(x)-1
        h = x(i+1)-x(i-1);
        df(i) = (y(i+1)-y(i-1))/h;
    end
end

%---------------------------Finds the Lagrange interpolation function
function [f,a] = Lagrange(xp,yp)
    p = length(xp);
    f = @(x)(0);
    for i = 1:p
        g = @(x)(1);
        for j = 1:p
            if i ~= j
                g = @(x)(g(x).*(x-xp(j))/(xp(i)-xp(j)));
            end
        end
        f = @(x)(f(x) + yp(i)*g(x));
    end
    xs = min(xp):0.01:max(xp);
    a = polyfit(xs,f(xs),p-1);
end

%---------------------------Finds the approximate slope at a point when given a function
function [DF] = Continuous_Diff(xp,f,h)
% Forward Difference
    fdf = (f(xp+h)-f(xp))./h;
% Backward Difference
    bdf = (f(xp-h)-f(xp))./-h;
% Centered Difference
    cdf = (f(xp+h)-f(xp-h))./(2*h);
% Table
    Data = [fdf; bdf; cdf];
    x1df = Data(:,1);
    x2df = Data(:,2);
    Cols = ["x_2";"x_2p5"];
    Rows = ["Forward","Backward","Centered"];
    DF = table(x1df,x2df,'VariableNames',Cols,'RowNames',Rows);
end

%% Problem 4 Functions
%---------------------------Gaussian Elimination Direct Solver
function [x] = Gaussian_Elim(A,b)
    n = length(b);
    % Forward Elimination
    for col = 1:n-1
        for row = col+1:n
            factor = A(row,col)/A(col,col);
            A(row,col) = 0;
            A(row,col+1:n) = A(row,col+1:n)-factor*A(col,col+1:n);
            b(row) = b(row)-factor*b(col);
        end
    end
    % Backward Substitution
    for i = n:-1:1
        temp = 0;
        for j = i+1:n
            temp = temp + A(i,j)*x(j);
        end
        x(i) = (b(i)-temp)/A(i,i);
    end
    x = x';
end