function output = btimes(a, b)
%% btimes Utility function shortening bsxfun(@times, a, b)
% Syntax:
%     output = btimes(a,b)
%     
% Comment:
%     Short for bsxfun(@times, a, b), removing the explicit use of @.

%   Created by: Rikard Vinge
%   $Revision: 1.0$  $Date: 2015/03/12 16:30:00$

output = bsxfun(@times, a, b);
end