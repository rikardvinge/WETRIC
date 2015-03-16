function output = bminus(a, b)
%% bminus Utility function shortening bsxfun(@minus, a, b)
% Syntax:
%     output = bminus(a,b)
%     
% Comment:
%     Short for bsxfun(@minus, a, b), removing the explicit use of @.

%   Created by: Rikard Vinge
%   $Revision: 1.0$  $Date: 2015/03/12 16:30:00$

output = bsxfun(@minus, a, b);
end