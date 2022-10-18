% Jacobi Method of Solving Ax=b for x by Kylie Hoyt
% An Iteratvie Method

close all
clear all

% Variable Definition
maxiter = 50;
tol = 1*10^-8;
err = 1.0;
iter = 0;
sor = 0;
lambda = 1;
A = [5 1 0; 2 9 3; 4 2 1];
b = [1 1 1]';
n = length(A);
xnew = ones(length(A),1);
xold = xnew;
A\b
% Convert Ax=b to xnew=d-C*xold where d=b/pivots and C=A/pivots with
% pivots set to 0
for i = 1:n
    b(i) = b(i)/A(i,i);
    A(i,:) = A(i,:)/A(i,i);
    A(i,i) = 0;
end
% Solve xnew=d-C*xold until err<tol using the xold vector for all xnew(i)
% calculations (vectorized)
while err>tol && iter<maxiter
   xnew = b-A*xnew;
   if sor == 1
       xnew = lambda*xnew+(1-lambda)*xold;
   end
   err = norm(xnew-xold,Inf)/norm(xnew,Inf);
   xold = xnew;
   iter = iter+1;
end
xnew