function output = brdivide(a, b)
%% brdivide Utility function shortening bsxfun(@rdivide, a, b)
% Syntax:
%     output = brdivide(a,b)
%     
% Comment:
%     Short for bsxfun(@rdivide, a, b), removing the explicit use of @.

%   Created by: Rikard Vinge
%   $Revision: 1.0$  $Date: 2015/03/12 16:30:00$

output = bsxfun(@rdivide, a, b);
end