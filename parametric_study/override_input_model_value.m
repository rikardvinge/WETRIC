function modelStructCell = override_input_model_value(inputModel, overrideParameter)
%% OVERRIDE_INPUT_MODEL_VALUE Plot utility for parametric results
%
% Syntax:
%     modelStructCell = OVERRIDE_INPUT_MODEL_VALUE(inputModel, overrideParameter)
%
% Input:
%     inputModel           - [cell OR string] name of input model(s)
%     overrideParameter    - [cell of struct OR struct] parameters to override
%
% Output:
%     modelStructCell      - [cell] cell array of ModelStruct data with
%                            parameterValues overrided according to
%                            overrideParameter input.
% Comment:
%     Utility function for simpler override of parameter values when computing
%     parametric studies.
%
% Example usage:
%     modelStructCell = override_input_model_value('analyt_collapsed_2coil', struct('N1',25))
%
% See also compute_parametric_study

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2015-03-16 14:00:00$

%% Check input:
if ischar(inputModel)
   inputModel = {inputModel};
end
if isstruct(overrideParameter)
   overrideParameter = {overrideParameter};
end
% Check for equal number of inputs:
modelNum = length(inputModel);
if length(overrideParameter)~=modelNum
   error('The number of inputModels and overrideParameter structures are not equal.')
end
%% Override parameterValues and generate ModelStruct input for parametric study:
modelStructCell           = cell(1,modelNum);
for modelIdx = 1:modelNum
   %% Load parameter structure in model:
   ModelStruct = load(inputModel{modelIdx});
   
   %% Check if single function call or multiple functions
   if ~iscell(ModelStruct.functionName)
      % Scalar case -> convert model struct to multiple array case
      ModelStruct.functionName    = {ModelStruct.functionName};
      ModelStruct.inputParameter  = {{},ModelStruct.inputParameter};
      ModelStruct.inputValue      = {{},ModelStruct.inputValue};
      ModelStruct.outputParameter = {ModelStruct.outputParameter};
   end
   
   % Check number of subfunction in the model:
   subFuncNum = length(ModelStruct.functionName);
   
   % Create functionSupportArrayInput field if it does not exist:
   if ~isfield(ModelStruct,'functionSupportArrayInput')
      ModelStruct.functionSupportArrayInput = false(1,subFuncNum);
   end
   
   % Override default values according to overrideParameter struct:
   if ~isempty(overrideParameter{modelIdx})
      overrideParName = fieldnames(overrideParameter{modelIdx});
      for overrideIdx = 1:length(overrideParName)
         flagInvalidParameter = true;
         for subFuncIdx = 1:subFuncNum+1
            % Find index of input parameter:
            tmpParIdx = find(strcmp(overrideParName{overrideIdx},...
               ModelStruct.inputParameter{subFuncIdx}));
            if ~isempty(tmpParIdx)
               flagInvalidParameter = false;
            else
               continue
            end
            % Override value according to field name:
            ModelStruct.inputValue{subFuncIdx}{tmpParIdx} = ...
               overrideParameter{modelIdx}.(overrideParName{overrideIdx});
         end
         if flagInvalidParameter
               error('Parameter "%s" is not a valid parameter for model "%s"',...
                  overrideParName{overrideIdx},inputModel{modelIdx})
               % Note, uncesessary control in a proper GUI
         end
      end
   end
   % Note, override structure is placeholder for a GUI where each value is set
   % in a separate box for to each input.
      
   % Save ModelStruct for these input models in a cell array with correct scalar
   % values:
   modelStructCell{modelIdx} = ModelStruct;
end