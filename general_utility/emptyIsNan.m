function out = emptyIsNan(inp)
%% emptyIsNan
% Syntax:
%  out = emptyIsNan(inp)
%
% Comments:
%  Converts empty elements in inp to nan values
tf       = isempty(inp);
out      = inp;
out(tf)  = nan;
