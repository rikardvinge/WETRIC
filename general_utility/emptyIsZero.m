function out = emptyIsZero(inp)
%% emptyIsZero
% Syntax:
%  out = emptyIsZero(inp)
%
% Comments:
%  Converts empty elements in inp to nan values
tf       = isempty(inp);
out      = inp;
out(tf)  = 0;
