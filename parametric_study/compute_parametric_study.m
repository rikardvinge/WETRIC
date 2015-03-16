function [ParametricResult, ParametricDimHelper] = compute_parametric_study(modelStructCell, ParametricStudy)
%% COMPUTE_PARAMETRIC_STUDY Computes a parametric study for the supplied model
%
% Syntax:
%     [ParametricResult, ParametricDimHelper] = COMPUTE_PARAMETRIC_STUDY(modelStructCell, ParametricStudy)
%
% Input:
%     modelStructCell      - [cell array] of ModelStruct data
%     ParametricStudy      - [structure] specifying the parametric study
%
% Output:
%     ParametricResult     - [structure] with fieldnames with input and output
%                            paramterers in the model(s) with corresponding
%                            values with dimensions depending on the parametric
%                            study input.
%     ParametricDimHelper  - [structure] with dimension size of parametric study
%                            and helpfull select/expand index arrays.
% Comment:
%
% Example usage:
%
% See also merge_struct, emptyIsZero

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2015-03-12 09:00:00$

% Control input:
modelNum                   = length(modelStructCell);

% Parse ParametricStudy structure and generate ND-matrix array of input
parametricParameter  = fieldnames(ParametricStudy);

% Read parametric study structure and generate ND-matrices
parametricParameterNum        = length(parametricParameter);
parametricParameterDim        = cellfun(@(x) ParametricStudy.(x){1},parametricParameter);
[parametricParameterDim,sortIdx]   = sort(parametricParameterDim);
[parametricParameterDimUniq,ia,ic] = unique(parametricParameterDim);
% parametricParameterDimNum     = max(parametricParameterDim);

% Sort parameter such that their dimensions are in increasing order:
parametricParameter           = parametricParameter(sortIdx);

% Update parameter dimensions to the lowest applicable dimension range:
parametricParameterDimRange   = 1:length(parametricParameterDimUniq);
for paramIdx = 1:parametricParameterNum
   ParametricStudy.(parametricParameter{paramIdx}){1} = ...
      parametricParameterDimRange(ic(paramIdx));
end
% Update
parametricParameterDim        = cellfun(@(str) ParametricStudy.(str){1},parametricParameter);
parametricParameterDimNum     = max(parametricParameterDim);
[parametricParameterDimUniq,ia,ic] = unique(parametricParameterDim);

% Generate ndgrid over parametric variables:
parametricParameterLinear     = cellfun(@(str) ParametricStudy.(str){2},parametricParameter,'un',0);
parametricParameterLinear     = cellfun(@(lval,ldim) shiftdim(lval(:),1-ldim),...
   parametricParameterLinear,num2cell(parametricParameterDim),'un',0);
parametricParameterSize       = cellfun(@(lval) length(lval),parametricParameterLinear);
% [sizeUniq,ia2,ic2] = unique(sizePar,'stable');
parametricParameterSizeUniq   = parametricParameterSize(ia);
parametricParameterSizeAll    = ones(1,parametricParameterDimNum);
parametricParameterSizeAll(parametricParameterDimUniq) = parametricParameterSizeUniq;

parametricParameterNumberOfComputations = prod(parametricParameterSizeUniq);

% Check memory size of a single nd-grid:
parametricParameterMemoryUsage = parametricParameterNumberOfComputations*8/1024^3;
if parametricParameterMemoryUsage*2 > 2
   error('The parametric study memory size is on the order of 2 GB.')
end

parametricParameterSelectAll  = arrayfun(@(dimSize) 1:dimSize,parametricParameterSizeAll,'un',0);
parametricParameterExpandAll  = arrayfun(@(dimSize) ones(1,dimSize),parametricParameterSizeAll,'un',0);

% Start allocation of results structure:
parametricParameterLinear  = cat(1,parametricParameter(:).',parametricParameterLinear(:).');
ParametricResult           = struct(parametricParameterLinear{:});
ParametricDimHelper        = struct('sizeAll',parametricParameterSizeAll);
ParametricDimHelper.selectAll = parametricParameterSelectAll;
ParametricDimHelper.expandAll = parametricParameterExpandAll;

% Clear some large variables
clear parametricParameterLinear
%% Loop over input models and compute output data for the parameteric study
for modelIdx = 1:modelNum
   
   %% Select current model structure:
   ModelStructExt = modelStructCell{modelIdx};
   
   % Check number of subfunction in the model:
   subFuncNum = length(ModelStructExt.functionName);
   
   %% Loop over each sub function in the model
   for subFuncIdx = 1:subFuncNum
      
      % Generate simple 'scalar' ModelStruct from ModelStructExt:
      ModelStruct = ModelStructExt;
      ModelStruct.functionName   = ModelStruct.functionName{subFuncIdx};
      ModelStruct.inputParameter = cat(2,ModelStruct.inputParameter{[1,subFuncIdx+1]});
      ModelStruct.inputValue     = cat(2,ModelStruct.inputValue{[1,subFuncIdx+1]});
      ModelStruct.outputParameter = ModelStructExt.outputParameter{subFuncIdx};
            
      %% Find and replace inputParameter to model already in ParametricResults
      % Note, this automatically replaces input values according to the parametric
      % study.
      
      % Update current parameters in the parametric results structure:
      parametricResultParameter = fieldnames(ParametricResult);
      
      parametricResultIdxInModel = cellfun(@(tmpParParam) emptyIsZero( ...
         find(strcmp(tmpParParam,ModelStruct.inputParameter)) ),...
         parametricResultParameter);
      parametricResultIsInModel = parametricResultIdxInModel ~= 0;
      parametricResultIdxInModel = parametricResultIdxInModel(parametricResultIsInModel);
      
      % Update for parametric parameters:
      for resultParameterIdx = 1:length(parametricResultIdxInModel)
         % Update model inputValues according to parametric parameter grid:
         ModelStruct.inputValue{parametricResultIdxInModel(resultParameterIdx)} = ...
            ParametricResult.(...
            ModelStruct.inputParameter{parametricResultIdxInModel(resultParameterIdx)});
      end
      
      %% Find and replace linked parameter values:
      
      % Step through all inputValue and check for strings with @signs:
      for inputParamIdx = 1:length(ModelStruct.inputParameter)
         if ~ischar(ModelStruct.inputValue{inputParamIdx})
            continue
         end
         if isempty(strfind(ModelStruct.inputValue{inputParamIdx},'@'))
            continue
         end
         tmpStr = ModelStruct.inputValue{inputParamIdx};
         tmpStr = strrep(tmpStr,'@','ParametricResult.');
         % Update model inputValues according to parametric parameter grid:
         tmpData = eval(tmpStr);
         ModelStruct.inputValue{inputParamIdx} = tmpData;
      end
      
      %% Place all model inputParameters into ParametricResults
      tmpInputDoesAlreadyExist = cellfun(@(inputParameter) ...
         any(strcmp(inputParameter,fieldnames(ParametricResult))), ...
         ModelStruct.inputParameter);
      
      tmpInputParameter = cat(1,ModelStruct.inputParameter(~tmpInputDoesAlreadyExist),...
         ModelStruct.inputValue(~tmpInputDoesAlreadyExist));
      tmpInputParameter = struct(tmpInputParameter{:});
      
      ParametricResult = merge_struct(ParametricResult,tmpInputParameter);
      
      %% Determine which value dimensions are non-singleton in the function call
      % Expand to minimum ND-size necessary for all input values
      
      % Find which inputValues that are scalar:
      tmpInpValueIsScalar = cellfun(@(inpParValue) isscalar(inpParValue),ModelStruct.inputValue);
      
      % Determine active dimension from the non-scalar inputValues:
      tmpInpValueSize      = cellfun(@(inpParValue) size(inpParValue),ModelStruct.inputValue,'un',0);
      tmpInpValueActiveDim = cellfun(@(inpParSize) emptyIsZero(find(inpParSize > 1)),tmpInpValueSize,'un',0);
      
      tmpInpValueUsedDim   = unique(cat(2,tmpInpValueActiveDim{:}));
      tmpInpValueUsedDim(tmpInpValueUsedDim == 0) = [];
      
      modelDimSelect       = num2cell(ones(1,parametricParameterDimNum));
      modelDimSelect(tmpInpValueUsedDim) = ...
         parametricParameterSelectAll(tmpInpValueUsedDim);
      
      % Expand input values according to the non-singleton dimensions:
      tmpInpValueNonScalarIdx = find(~tmpInpValueIsScalar);
      
      for inputIdx = 1:length(tmpInpValueNonScalarIdx)
         % Check active dimension in this input field:
         tmpParActiveDim = tmpInpValueActiveDim{tmpInpValueNonScalarIdx(inputIdx)};
         % Expand values to full modelDimSelect range:
         tmpExpandDim = setdiff(tmpInpValueUsedDim,tmpParActiveDim);
         tmpModelDimSelectExpand = modelDimSelect;
         tmpModelDimSelectExpand(tmpExpandDim) = ...
            parametricParameterExpandAll(tmpExpandDim);
         ModelStruct.inputValue{tmpInpValueNonScalarIdx(inputIdx)} = ...
            ModelStruct.inputValue{tmpInpValueNonScalarIdx(inputIdx)}( ...
            tmpModelDimSelectExpand{:} );
      end
      
      %% Compute outputValues according to inputValues using model function      
      
      % Preallocate cell array:
      modelOutputValue = cell(1,length(ModelStruct.outputParameter));
         
      
      if ModelStruct.functionSupportArrayInput(subFuncIdx)
         % Compute using ND-array data in call to function:
         [ modelOutputValue{:} ] = feval(ModelStruct.functionName, ...
            ModelStruct.inputValue{:});
      else
         % Compute using for loop over scalar data in call to function:
         
         % Arrayfun version [scalar input at a time - large overhead?!]:          
         % We can use arrayfun if we expand all the inputs:
         modelExpandDim = num2cell(ones(1,parametricParameterDimNum));
         modelExpandDim(tmpInpValueUsedDim) = ...
            parametricParameterExpandAll(tmpInpValueUsedDim);
         % Expand input values according to the non-singleton dimensions:
         tmpInpValueScalarIdx = find(tmpInpValueIsScalar);

         for inputIdx = 1:length(tmpInpValueScalarIdx)
            % Expand the scalar input:
            ModelStruct.inputValue{tmpInpValueScalarIdx(inputIdx)} = ...
               ModelStruct.inputValue{tmpInpValueScalarIdx(inputIdx)}( ...
               modelExpandDim{:} );
         end
         
         % Arrayfun over all [extended] input:
         func = str2func(ModelStruct.functionName);
         [ modelOutputValue{:} ] = arrayfun(func,...
            ModelStruct.inputValue{:});
         % Note, there is a significant slow down if we call a function with
         % many small sub-functions.
         
%          % Alternate version [explicit for loop over all inputs - similar]:
%          computeNum = prod(cellfun(@(select) length(select),modelExpandDim));
%          
%          modelOutputValue = repmat({zeros(size(ModelStruct.inputValue{1}))},1,length(ModelStruct.outputParameter));
%          numOutputParam  = length(ModelStruct.outputParameter);
%          tmpScalarOutput = cell(1,numOutputParam);
%          for computeIdx = 1:computeNum
%             
%             tmpScalarInput = cellfun(@(data) data(computeIdx),...
%                ModelStruct.inputValue,'un',0);
%             
%             [tmpScalarOutput{:}] = feval(ModelStruct.functionName,...
%                tmpScalarInput{:});
%             
%             for outputIdx = 1:numOutputParam
%                modelOutputValue{outputIdx}(computeIdx) = ...
%                   tmpScalarOutput{outputIdx};
%             end            
%          end
         
         
      end
      
      % Find which outputValues are scalar:
      tmpOutValueIsScalar     = cellfun(@(outParValue) ...
         isscalar(outParValue),modelOutputValue);
      
      % Check the non-scalar output values for their dependence:
      tmpOutValueNonScalarIdx = find(~tmpOutValueIsScalar);
      tmpAllDim               = 1:parametricParameterDimNum;
      
      for outputIdx = 1:length(tmpOutValueNonScalarIdx)         
         % Check absolute sum difference over each dimension:
         tmpSumAbsDiff = arrayfun(@(dim) ...
            sum(reshape(abs(diff(modelOutputValue{...
            tmpOutValueNonScalarIdx(outputIdx)},1,dim)),[],1)),tmpAllDim);
         
         % Check for which dimensions we have no difference
         tmpIsNotDependent = tmpSumAbsDiff == 0; % < 10*(eps);
         % [TODO, relative check?]: Find mean value of non-scalar outputs:
         %    tmpOutMeanValue     = cellfun(@(outValue) mean(outValue(:)), )
         
         % Store only for dimension where there are any difference:
         tmpModelDimStore = modelDimSelect;
         tmpModelDimStore(tmpIsNotDependent) = ...
            num2cell(ones(1,nnz(tmpIsNotDependent)));
         
         modelOutputValue{tmpOutValueNonScalarIdx(outputIdx)} = ...
            modelOutputValue{tmpOutValueNonScalarIdx(outputIdx)}( ...
            tmpModelDimStore{:} );
      end
      
      % Place all model outputParameters into ParametricResults
      modelOutputValue = cat(1,ModelStruct.outputParameter, modelOutputValue);
      modelOutputValue = struct(modelOutputValue{:});
      
      ParametricResult = merge_struct(ParametricResult,modelOutputValue);
      clear tmpModelOutput      
   end
end