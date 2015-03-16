function output = bplus(a, b)
%% bplus Utility function shortening bsxfun(@plus, a, b)
% Syntax:
%     output = btimes(a,b)
%     
% Comment:
%     Short for bsxfun(@plus, a, b), removing the explicit use of @.

%   Created by: Rikard Vinge
%   $Revision: 1.0$  $Date: 2015/03/12 16:30:00$

output = bsxfun(@plus, a, b);
end