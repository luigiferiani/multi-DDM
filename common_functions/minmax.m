function [ extremes ] = minmax( A )
%minmax Gives minimum and maximum value of a matrix

if ~isnumeric(A)
    error('A has to be numeric');
end

m = min(A(:));
M = max(A(:));

extremes = [m, M];


end

