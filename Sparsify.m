function [D,loc] = sparsify(A)
n = length(A);
loc = [];
for i = 1:n
    if A(i,1) ~= 0 
        loc(length(loc)+1) = -i+1;
    end
end
for j = 2:n
    if A(1,j) ~= 0 
        loc(length(loc)+1) = j-1;
    end
end
loc = sort(loc);
D = zeros(n,length(loc));
Alarge = [zeros(n,abs(min(loc))),A,zeros(n,max(loc))];
shiftedlocs = loc + abs(min(loc))+1;
for i = 1:n
    D(i,:) = Alarge(i,i+shiftedlocs-1);
end
end
