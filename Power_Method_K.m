% Power Method of Finding Eigenvectors and Eigenvalues by Kylie Hoyt

maxiter = 100;
tol = 0.01;
xold = [1 1]';
xnew = xold;
err = 1.0;
iter = 0;
A = [0 -1; 2 3];

while err>tol && iter<maxiter
    iter = iter+1;
    xnew = A*xnew;
    xnew = xnew/norm(xnew);
    err = norm(xnew-xold)/norm(xnew);
    lambda = xnew'*A*xnew;
    xold = xnew;
end
lambda
xnew