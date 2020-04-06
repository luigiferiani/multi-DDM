function [] = setsemilogx( axis_handle )
%setsemilogx Set the axis in semilogarithmic scale

if nargin < 1 || isempty(axis_handle)
    axis_handle = gca;
end

set(axis_handle,'XScale','log','YScale','lin');


end

