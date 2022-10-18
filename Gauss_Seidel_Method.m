% Gauss Seidel Method of Solving Ax=b for x by Kylie Hoyt
% An Iterative Method

% Variable Definition
maxiter = 50;
tol = 1*10^-5;
err = 1.0;
iter = 0;
sor = 0; % Variant Flag
lambda = 1; % Variant
A = [5 1 0; 2 9 3; 4 2 1];
b = [1 1 1]';
n = length(b);

% Convert Ax=b to xnew=d-C*xold where d=b/pivots and C=A/pivots with
% pivots set to 0 (non-vectorized but converges quicker)
for i = 1:n
    b(i) = b(i)/A(i,i);
    A(i,:) = A(i,:)/A(i,i);
    A(i,i) = 0;
end
xnew = zeros(size(b));
xold = xnew;
% Solve xnew=d-C*xold until err<tol using the most recent xnew(i)
while err>tol && iter<maxiter
   iter = iter+1;
   for i = 1:n
       xnew(i) = b(i)-A(i,:)*xnew;
   end 
   if sor == 1
       xnew = lambda*xnew+(1-lambda)*xold;
   end
   err = norm(xnew-xold,Inf)/norm(xnew,Inf);
   xold = xnew;
end
iter
xnew