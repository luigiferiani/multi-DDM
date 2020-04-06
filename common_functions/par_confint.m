function [ ci ] = par_confint( fit_out, par_name, level )
%par_confint return confidence interval for the parameter in par_name
%   Detailed explanation goes here

if ~ischar(par_name)
    error('need chars in input');
end
if nargin < 3 || isempty(level)
    level = 0.95;
end

% list the coefficient names
parnames = coeffnames(fit_out);

% find the one that matches the input
parind = find(strcmp(parnames,par_name));

if isempty(parind)
    error('parameter not found, check your input.');
end

% get confindence interval
cimat = confint(fit_out, level);

ci = cimat(:,parind);

end

