%% Lagrange Interpolation

function [f] = Lagrange(x,y,n)
    f = 0;
    for i = 1:n
        g = 1;
        for j = 1:n
            if i ~= j
                g = @(x)(g*(x-x(j))/(x(i)-x(j)));
            end
        end
        f = @(x)(f + y(i)*g);
    end
end