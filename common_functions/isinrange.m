function [ idx_inrange ] = isinrange( x, minvalue, maxvalue )
%ISINRANGE returns true if x is within minvalue and maxvalue. 

if (nargin <3 || isempty(maxvalue)) && numel(minvalue) == 2
    maxvalue = max(minvalue);
    minvalue = min(minvalue);
elseif (nargin == 3 && ~isempty(maxvalue)) && (numel(minvalue) > 1 || numel(maxvalue) > 1)
    minvalue = min(minvalue(:));
    maxvalue = max(maxvalue(:));
end %if

mM = minmax([minvalue, maxvalue]);
minvalue = mM(1);
maxvalue = mM(2);

idx_inrange = (x >= minvalue) & (x <= maxvalue);


end

