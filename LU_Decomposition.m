clear all

close all

% LU Decomposition: A direct method to solve Ax=b for many b's

% Initialize Matrices
A = [2 1 -1;-3 -1 2;-2 1 2];
b = [8 -11 -3;7 4 9;3 0 6]';
n = length(A);

n = 4; % number of elements
e = ones(n,1);
b = zeros(n,1);
b(1,1) = 1;
b(n,1) = -1;
sp_st = [e -2*e e];
A = spdiags(sp_st,[-2 0 1],n,n);

% Forward Elimination: A = LU
for col = 1:n-1
    for row = col+1:n
        factor = A(row,col)/A(col,col);
        A(row,col) = 0;
        A(row,col+1:n) = A(row,col+1:n)-factor*A(col,col+1:n);
        A(row,col) = factor;
    end
end
full(A)
for k = 1:length(b)
    % Forward Substitution: Solve Ly = b for y
    for i = 1:n
        temp = 0;
        for j = 1:i-1
            temp = temp + A(i,j)*y(j);
        end
        y(i) = b(i,k)-temp;
    end

    % Backward Substitution: Solve Ux = y for x
    for i = n:-1:1
        temp = 0;
        for j = i+1:n
            temp = temp + A(i,j)*x(j);
        end
        x(i) = (y(i)-temp)/A(i,i);
    end
    x'
end
