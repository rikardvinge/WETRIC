function str = cellstr2str(cellstr, separationStr, numericConversionStr)
%% CELLSTR2STR Convert a cellstring to a single string with a separation string
%
% Syntax:
%     str = cellstr2str(cellstr, separationStr)
%     str = cellstr2str(numericArray, separationStr, numericConversionStr)
%
% Comment:
%     Utility function for writing a cellstr as a single string.
%     Can also convert a numeric array to a string.
%     Main purpose is add separation charachter between the parts of the
%     cellstr.

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2012/?/? 00:00:00$
%   $Revision: 2.0$  $Date: 2014/10/06 13:00:00$
%     Complete overhaul of function

% Default separationStr
if (nargin <= 1)
   separationStr = ' ';
end
% Default numericConversionStr
if nargin <= 2
   numericConversionStr = '%g';
end

% Numeric input array:
if ~iscell(cellstr)
   % Convert numeric array to cellstr:
   cellstr = cellfun(@(innum) num2str(innum, numericConversionStr),...
      num2cell(cellstr),'un',0);
end

% Concatenate into output:
cellstr        = cellstr(:).';
len_cs         = length(cellstr);
cellsty        = cat(2, repmat( { separationStr }, [1, len_cs-1]), {''});
cellstr_merge  = cat(1, cellstr, cellsty );
str            = cat(2, cellstr_merge{:} ) ;
