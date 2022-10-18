% Deflation Method of Finding Eigenvectors and Eigenvalues by Kylie Hoyt

maxiter = 100;
tol = 0.00001;
err = 1.0;
iter = 0;
A = [7 4; -6 -3];

for i = 1:size(A)
    xold = ones(length(A),1);
    xnew = xold;
    lambda = 1;
    while err>tol && iter<maxiter
        iter = iter+1;
        xnew = A*xnew;
        xnew = xnew/norm(xnew);
        err = norm(A*xnew-lambda*xnew);
        lambda = xnew'*A*xnew;
        xold = xnew;
    end
    lambda
    xnew
    A = A-lambda*(xnew)*(xnew');
    err = 1.0;
    iter = 0;
end