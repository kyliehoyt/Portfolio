

close all

% Gaussian Elimination: A direct method to solve Ax=b
% Initialize Matrices
A = [5 1 0; 2 9 3; 4 2 1];
b = [1 1 1]';
n = length(A);

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
x'
