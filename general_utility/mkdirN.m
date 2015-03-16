function mkdirN(folder)
%% mkdirN Utility function creating a directory unless it already exists
% Syntax:
%     mkdirN(folder)
%     
% Comment:
%     Creates a folder if it does not already exist!

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2012/?/? 00:00:00$
%   $Revision: 2.0$  $Date: 2014/10/06 13:00:00$
%     Cleaned function script and comments

if ~exist(folder,'dir')
	mkdir(folder)
end
