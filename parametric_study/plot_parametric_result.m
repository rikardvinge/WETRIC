function plot_parametric_result(fileName, pathName)
%% PLOT_PARAMETRIC_RESULT Plot utility for parametric results
%
% Syntax:
%     PLOT_PARAMETRIC_RESULT()
%     PLOT_PARAMETRIC_RESULT(fileName, pathName)
%
% Input:
%     fileName             - [string] name of parametric result file
%     pathName             - [string] path to parametric result file
%
% Output:
%     plot utility GUI figure for regular line plots
%
% Comment:
%
% Example usage:
%
% See also emptyIsZero

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2015-03-17 15:00:00$

%% Load input file:

% Use uigetfile if not input
if nargin == 0
   [fileName,pathName,filterIndex] = ...
      uigetfile('.mat','Select parametric result file','..\WETRICdata');
end

% Load results:
resultData = load([pathName fileName]);

%% Parse ParametricResult parameters
parameterStr      = fieldnames(resultData.ParametricResult);
ParameterSize     = structfun(@(data) size(data), resultData.ParametricResult,'un',0);
parameterDimNum   = max(structfun(@(data) emptyIsZero(find(data~=1,1,'last')), ParameterSize));
parameterIsScalar = structfun(@(data) all(data==1), ParameterSize);
parameterIsVector = structfun(@(data) nnz(data~=1)==1, ParameterSize);
parameterIsGrid   = ~(parameterIsScalar | parameterIsVector);

ParameterActiveDim = structfun(@(data) emptyIsZero(find(data~=1)), ParameterSize,'un',0);

parameterStrIfVector = parameterStr(parameterIsVector);
parameterActiveDimIfVector = cellfun(@(str) ParameterActiveDim.(str),parameterStrIfVector);

parameterStrMaxLen = max(cellfun(@(str) length(str), parameterStr));
tmpStrFormat       = ['%-' sprintf('%d',parameterStrMaxLen) 's = %s'];
precisionDigits    = 3;
% Generate text of selected scalar values
textScalarParam  = cellfun(@(param) sprintf(tmpStrFormat, param, ...
   num2sci(resultData.ParametricResult.(param), precisionDigits)), ...
   parameterStr(parameterIsScalar), 'un', 0);

% Generate text with size of grid parameters
textGridParam    = cellfun(@(param) sprintf(tmpStrFormat,param,...
   mat2str(ParameterSize.(param))),...
   parameterStr(parameterIsGrid),'un',0);

% Generate string array of the vectorial parameters
[parameterActiveDimIfVector,sortIdx] = sort(parameterActiveDimIfVector);
textVectParam    = parameterStr(parameterIsVector);
textVectParam     = textVectParam(sortIdx);
textVectParamDim = arrayfun(@(dim) sprintf('dim %d: ',dim),parameterActiveDimIfVector,'un',0);

textVectParamInformative = cellfun(@(strDim,strPar) cat(2,strDim,strPar),textVectParamDim,textVectParam,'un',0);

% Set up plot options
PlotOptionDef              = struct();
PlotOptionDef.xParam.desc  = 'x data';
PlotOptionDef.xParam.type  = 'popupmenu';
PlotOptionDef.xParam.value = textVectParamInformative;
PlotOptionDef.xParam.data  = textVectParam;


PlotOptionDef.xScale.desc  = 'x scaling';
PlotOptionDef.xScale.type  = 'text';
PlotOptionDef.xScale.value = 1;
PlotOptionDef.xScale.data  = 'num';

PlotOptionDef.xLabel.desc  = 'x label';
PlotOptionDef.xLabel.type  = 'text';
PlotOptionDef.xLabel.value = ''; %'x [-]';
PlotOptionDef.xLabel.data  = 'str';

PlotOptionDef.xLogScale.desc   = 'x logscale';
PlotOptionDef.xLogScale.type   = 'checkbox';
PlotOptionDef.xLogScale.value  = false;
PlotOptionDef.xLogScale.data   = 'num';

PlotOptionDef.yExpr.desc   = 'y expr';
PlotOptionDef.yExpr.type   = 'text';
PlotOptionDef.yExpr.value  = ''; %'real(@vRl.*conj(@iRl)./(@vC1.*conj(-@iRg)))';
PlotOptionDef.yExpr.data  = 'str';

PlotOptionDef.yLabel.desc  = 'y label';
PlotOptionDef.yLabel.type  = 'text';
PlotOptionDef.yLabel.value = ''; %'y [-]';
PlotOptionDef.yLabel.data  = 'str';

PlotOptionDef.yLogScale.desc   = 'y logscale';
PlotOptionDef.yLogScale.type   = 'checkbox';
PlotOptionDef.yLogScale.value  = false;
PlotOptionDef.yLogScale.data   = 'num';

PlotOptionDef.yMinMaxRange.desc   = 'data range';
PlotOptionDef.yMinMaxRange.type   = 'text';
PlotOptionDef.yMinMaxRange.value  = '[-inf, inf]'; %'[0, 1]';
PlotOptionDef.yMinMaxRange.data   = 'num';

PlotOptionDef.legend.desc   = 'Legend if>1';
PlotOptionDef.legend.type   = 'checkbox';
PlotOptionDef.legend.value  = true;
PlotOptionDef.legend.data   = 'num';

PlotOptionDef.title.desc   = 'title';
PlotOptionDef.title.type   = 'text';
PlotOptionDef.title.value  = '';
PlotOptionDef.title.data  = 'str';

PlotOptionDef.lw.desc   = 'line width';
PlotOptionDef.lw.type   = 'text';
PlotOptionDef.lw.value  = '1';
PlotOptionDef.lw.data  = 'num';

PlotOptionDef.lt.desc   = 'line style';
PlotOptionDef.lt.type   = 'text';
PlotOptionDef.lt.value  = '-';
PlotOptionDef.lt.data  = 'str';

plotOptionStr = fieldnames(PlotOptionDef);
plotOptionNum = length(plotOptionStr);

% UI-control options:
GuiOpt = struct();
GuiOpt.fontName               = 'times';
GuiOpt.textFontSize           = 11;
GuiOpt.editFontSize           = 10;
GuiOpt.popupFontSize          = 9;
GuiOpt.pushbuttonFontSize     = 12;

GuiOpt.pixelHeigthUIcontrol   = 20;
GuiOpt.pixelPaddingHeigth     = [8,   5]; % [bottom/top, between uicontrols]
GuiOpt.pixelPaddingWidth      = [6,   4]; % [bottom/top, between uicontrols]

GuiOpt.extentWidthOptionText    = 90;  % It is increased when necessary
GuiOpt.extentWidthOptionEdit    = 270; % The width of the edit field
GuiOpt.extentWidthTotal         = GuiOpt.extentWidthOptionText + GuiOpt.extentWidthOptionEdit +GuiOpt.pixelPaddingWidth(2);
GuiOpt.textAlignPad             = 2;
GuiOpt.numLines                 = 7;
GuiOpt.extentControlWindowText  = 80;
GuiOpt.extentControlWindowSlider = 200;
GuiOpt.extentControlWindowTotal = GuiOpt.extentControlWindowText + GuiOpt.extentControlWindowSlider + GuiOpt.pixelPaddingWidth(2);

% Generate a figure
hFig           = figure('Name','Plot Parametric Result',...
   'Toolbar','none','Menubar','none','NumberTitle','off');

% Print uicointrol components start from the bottom of the list:
tmpCurrentPosition  = [GuiOpt.pixelPaddingWidth(1), GuiOpt.pixelPaddingHeigth(1)];

% Create "Plot" button
tmpPosition   = [tmpCurrentPosition, ...
   GuiOpt.extentWidthOptionText, GuiOpt.pixelHeigthUIcontrol];
hPlotButton   = uicontrol( ...
   'Style',       'Pushbutton',...
   'String',      'Plot',...
   'FontName',    GuiOpt.fontName, ...
   'FontSize',    GuiOpt.pushbuttonFontSize, ...
   'Position',    tmpPosition);
% Check actual necessary size:
tmpExtent = get(hPlotButton,'extent');
% extentWidthUniversal   = max(extentWidthUniversal, ...
%    tmpExtent(3)+extentWidthPushbuttonInPadding);
% tmpNewSize  = [extentWidthUniversal, ...
%    tmpExtent(4) + extentHeigthPushbuttonInPadding];
% Update size of uicontrol:
% set(hPlotButton,'Position',[tmpCurrentPosition,tmpNewSize]);
% Update current position Y coordinate:
tmpCurrentPosition(2) = tmpCurrentPosition(2) + tmpExtent(4) + ...
   GuiOpt.pixelPaddingHeigth(2);

% Create "Plot Option" uicontrols
hOptions       = cell(plotOptionNum,2);
for plotOptionIdx = plotOptionNum:-1:1
   
   % Get current plot option:
   tmpOption      = PlotOptionDef.(plotOptionStr{plotOptionIdx});
   
   % Generate current plot option size:
   tmpPosition    = [tmpCurrentPosition, ...
      GuiOpt.extentWidthOptionText, GuiOpt.pixelHeigthUIcontrol];
   
   % Create "desc" label text uicontrol:
   tmpPosition(2) = tmpPosition(2) - GuiOpt.textAlignPad;
   hOptions{plotOptionIdx,1} = uicontrol( ...
      'Style',       'text',...
      'String',      tmpOption.desc,...
      'FontName',    GuiOpt.fontName, ...
      'FontSize',    GuiOpt.textFontSize, ...
      'Position',    tmpPosition, ...
      'HorizontalAlignment', 'left');
   
   % Check actual necessary size:
   tmpExtent = get(hOptions{plotOptionIdx,1},'extent');
   %    GuiOpt.extentWidthOptionText   = ...
   %       max(GuiOpt.extentWidthOptionText, tmpExtent(3) );
   %    GuiOpt.pixelHeigthUIcontrol   = ...
   %       max(GuiOpt.pixelHeigthUIcontrol, tmpExtent(4) );
   tmpNewSize  = [GuiOpt.extentWidthOptionText, GuiOpt.pixelHeigthUIcontrol];
   % Update size of uicontrol:
   tmpPosition = [tmpCurrentPosition,tmpNewSize];
   tmpPosition(2) = tmpPosition(2) - GuiOpt.textAlignPad;
   set(hOptions{plotOptionIdx,1},'Position',tmpPosition);
   
   % Generate current plot option size:
   tmpCurrentPosition    = tmpCurrentPosition + ...
      [GuiOpt.pixelPaddingWidth(2)+GuiOpt.extentWidthOptionText,0];
   
   % Generate current plot option size:
   tmpPosition    = [tmpCurrentPosition, ...
      GuiOpt.extentWidthOptionEdit, GuiOpt.pixelHeigthUIcontrol];
   
   % Check plot option type:
   if strcmp(tmpOption.type,'text')
      
      % Create uicontrol:
      hOptions{plotOptionIdx,2} = uicontrol( ...
         'Style',       'edit',...
         'String',      num2str(tmpOption.value),...
         'FontName',    GuiOpt.fontName, ...
         'FontSize',    GuiOpt.editFontSize, ...
         'Position',    tmpPosition, ...
         'HorizontalAlignment', 'right',...
         'UserData',    {plotOptionStr{plotOptionIdx}, tmpOption.data });
      
   elseif strcmp(tmpOption.type,'popupmenu')
      
      % Create uicontrol:
      hOptions{plotOptionIdx,2} = uicontrol( ...
         'Style',       'popupmenu', ...
         'String',      tmpOption.value,...
         'FontName',    GuiOpt.fontName, ...
         'FontSize',    GuiOpt.popupFontSize, ...
         'Position',    tmpPosition, ...
         'HorizontalAlignment', 'right',...
         'UserData',    {plotOptionStr{plotOptionIdx}, tmpOption.data });
      
   elseif strcmp(tmpOption.type,'checkbox')
      
      % Create uicontrol:
      hOptions{plotOptionIdx,2} = uicontrol( ...
         'Style',       'checkbox',...
         'String',      'Enable',...
         'Value',       tmpOption.value,...
         'FontName',    GuiOpt.fontName, ...
         'FontSize',    GuiOpt.popupFontSize, ...
         'Position',    tmpPosition, ...
         'HorizontalAlignment', 'right',...
         'UserData',    {plotOptionStr{plotOptionIdx}, tmpOption.data });
      
   end
   
   % Update current position to next line
   tmpCurrentPosition(1) = GuiOpt.pixelPaddingWidth(1);
   tmpCurrentPosition(2) = tmpCurrentPosition(2) + tmpNewSize(2) + ...
      GuiOpt.pixelPaddingHeigth(2);
   
end

% Generate current plot option size:
tmpPosition    = [tmpCurrentPosition-[0,GuiOpt.pixelPaddingHeigth(2)], ...
   GuiOpt.extentWidthTotal, GuiOpt.pixelHeigthUIcontrol];
uicontrol( ...
   'Style',       'text',...
   'String',      'Plot options:',...
   'FontName',    GuiOpt.fontName, ...
   'FontSize',    GuiOpt.pushbuttonFontSize, ...
   'FontWeight',  'bold', ...
   'Position',    tmpPosition, ...
   'HorizontalAlignment', 'left');
tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
   GuiOpt.pixelHeigthUIcontrol + ...
   GuiOpt.pixelPaddingHeigth(2);

% Create "grid" parameter controls

% Generate current plot option size:
tmpPosition    = [tmpCurrentPosition, ...
   GuiOpt.extentWidthTotal, GuiOpt.numLines*GuiOpt.pixelHeigthUIcontrol];

% tmpStrCell = cellfun(@(str) sprintf('%s\n',str),textGridParam,'un',0);
str = sprintf(cellstr2str(textGridParam, '\n'));
uicontrol( ...
   'Style',       'edit',...
   'String',      str,...
   'Min',         0,...
   'Max',         2,...
   'Enable',      'inactive',...
   'FontName',    'fixedwidth', ...
   'FontSize',    GuiOpt.popupFontSize, ...
   'Position',    tmpPosition, ...
   'HorizontalAlignment', 'left');

tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
   GuiOpt.numLines*GuiOpt.pixelHeigthUIcontrol;

% Generate current plot option size:
tmpPosition    = [tmpCurrentPosition, ...
   GuiOpt.extentWidthTotal, GuiOpt.pixelHeigthUIcontrol];
uicontrol( ...
   'Style',       'text',...
   'String',      'Grid parameter sizes:',...
   'FontName',    GuiOpt.fontName, ...
   'FontSize',    GuiOpt.pushbuttonFontSize, ...
   'FontWeight',  'bold', ...
   'Position',    tmpPosition, ...
   'HorizontalAlignment', 'left');
tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
   GuiOpt.pixelHeigthUIcontrol + ...
   GuiOpt.pixelPaddingHeigth(2);

% Create "scalar" parameter controls

% Generate current plot option size:
tmpPosition    = [tmpCurrentPosition, ...
   GuiOpt.extentWidthTotal, GuiOpt.numLines*GuiOpt.pixelHeigthUIcontrol];

str = sprintf(cellstr2str(textScalarParam, '\n'));
uicontrol( ...
   'Style',       'edit',...
   'String',      str,...
   'Min',         0,...
   'Max',         2,...
   'Enable',      'inactive',...
   'FontName',    'fixedwidth', ...
   'FontSize',    GuiOpt.popupFontSize, ...
   'Position',    tmpPosition, ...
   'HorizontalAlignment', 'left');

tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
   GuiOpt.numLines*GuiOpt.pixelHeigthUIcontrol;

% Generate current plot option size:
tmpPosition    = [tmpCurrentPosition, ...
   GuiOpt.extentWidthTotal, GuiOpt.pixelHeigthUIcontrol];
uicontrol( ...
   'Style',       'text',...
   'String',      'Scalar parameter values:',...
   'FontName',    GuiOpt.fontName, ...
   'FontSize',    GuiOpt.pushbuttonFontSize, ...
   'FontWeight',  'bold', ...
   'Position',    tmpPosition, ...
   'HorizontalAlignment', 'left');
tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
   GuiOpt.pixelHeigthUIcontrol + ...
   GuiOpt.pixelPaddingHeigth(2);

% Update size of window
screenSize           = get(0,'ScreenSize');
figureSize           = [GuiOpt.extentWidthTotal + 2*GuiOpt.pixelPaddingWidth(1),...
   tmpCurrentPosition(2) + GuiOpt.pixelPaddingHeigth(1) - GuiOpt.pixelPaddingHeigth(2)];
figurePosition       = 0.5*screenSize([3,4]) - 0.5*figureSize;
set(hFig,'Position',[figurePosition,figureSize]);


% Set callback to plot button:
set(hPlotButton,'Callback',{@create_plot,hOptions(:,2),resultData.ParametricResult,GuiOpt});

drawnow

function create_plot(hObj,event,hOptions,ParametricResult,GuiOpt)
% Generate PlotOption struct from hOptions:
PlotOption = struct();
numOption = length(hOptions);
for optionIdx = 1:numOption
   
   % Get Option name and data type:
   optionStr      = hOptions{optionIdx}.UserData{1};
   optionDataType = hOptions{optionIdx}.UserData{2};
   
   % Get UI control style and values:
   optionUIstyle  = hOptions{optionIdx}.Style;
   optionUIvalue  = hOptions{optionIdx}.Value;
   optionUIstring = hOptions{optionIdx}.String;
   
   % Check plot option type:
   if strcmp(optionUIstyle,'edit')
      if strcmp(optionDataType,'num')
         optionValue = str2num(optionUIstring);
      elseif strcmp(optionDataType,'str')
         optionValue = optionUIstring;
      end
   elseif strcmp(optionUIstyle,'popupmenu')
      optionValue = optionDataType{optionUIvalue};
   elseif strcmp(optionUIstyle,'checkbox')
      if strcmp(optionDataType,'num')
         optionValue = optionUIvalue;
      end
   end
   
   PlotOption.(optionStr) = optionValue;
end

% Find all active dimensions for the parameters:
parameterStr      = fieldnames(ParametricResult);
ParameterSize     = structfun(@(data) size(data), ParametricResult,'un',0);
parameterDimNum   = max(structfun(@(data) emptyIsZero(find(data~=1,1,'last')), ParameterSize));
for parameterIdx = 1:length(parameterStr)
   tmpOnes = ones(1,parameterDimNum);
   tmpSize = ParameterSize.(parameterStr{parameterIdx});
   tmpOnes(1:length(tmpSize)) = tmpSize;
   ParameterSize.(parameterStr{parameterIdx}) = tmpOnes;
end
parameterIsVector = structfun(@(data) nnz(data~=1)==1, ParameterSize);
ParameterActiveDim = structfun(@(data) emptyIsZero(find(data~=1)), ParameterSize,'un',0);

parameterStrIfVector = parameterStr(parameterIsVector);
parameterActiveDimIfVector = cellfun(@(str) ParameterActiveDim.(str),parameterStrIfVector);

parameterSizeAll = struct2cell(ParameterSize);
parameterSizeAll = max(cat(1,parameterSizeAll{:}),[],1);

selectedIdxInEachDim = round(0.5*parameterSizeAll);

% Find X dimensions from ParametericStudy and PlotOption
xDim = find(size(ParametricResult.(PlotOption.xParam))>1);
% yDim = find(size(ParametricResult.(PlotOption.yParam))>1);

% Auto-generate x label, y label and title if empty:
if isempty(PlotOption.xLabel)
   PlotOption.xLabel = PlotOption.xParam;
end
if isempty(PlotOption.yLabel)
   PlotOption.yLabel = PlotOption.yExpr;
end
if isempty(PlotOption.title)
%    PlotOption.title = PlotOption.zExpr;
end


% Create selection array
parametricParameterExpandAll  = arrayfun(@(dimSize) ones(1,dimSize), parameterSizeAll,'un',0);
parametricParameterSelectAll  = arrayfun(@(dimSize) 1:dimSize, parameterSizeAll,'un',0);
selectIdxCellArray            = num2cell(selectedIdxInEachDim);
selectIdxCellArray(xDim)      = parametricParameterSelectAll(xDim);
% selectIdxCellArray(yDim)      = parametricParameterSelectAll(yDim);

% Select x data:
allDim = 1:length(parameterSizeAll);
permuteArray = unique(cat(2,xDim,allDim),'stable');
selectAxisIdxCellArray            = num2cell(ones(1,length(selectedIdxInEachDim)));
selectAxisIdxCellArray(xDim)      = parametricParameterSelectAll(xDim);
% selectAxisIdxCellArray(yDim)      = parametricParameterExpandAll(yDim);
xData = permute(ParametricResult.(PlotOption.xParam)(selectAxisIdxCellArray{:}),permuteArray);
% selectAxisIdxCellArray            = num2cell(ones(1,length(selectedIdxInEachDim)));
% selectAxisIdxCellArray(yDim)      = parametricParameterSelectAll(yDim);
% selectAxisIdxCellArray(xDim)      = parametricParameterExpandAll(xDim);
% yData = permute(ParametricResult.(PlotOption.yParam)(selectAxisIdxCellArray{:}),permuteArray);

% Evaluate zExpr:
tmpStr = PlotOption.yExpr;
% TODO? add functionalty of plotting multiple lines in the same plot!
% sepIdx = emptyIsZero( strfind(tmpStr,'|') );
% lineNum = sepIdx + 1;
% Split tmpStr into cell array at| sign (if any)
tmpStrEval = strsplit(tmpStr,'|');
lineNum = length(tmpStrEval);
yDataFullNDcell = cell(lineNum,1);
yDataLineName = cell(lineNum,1);
strMinCell = cell(lineNum,1);
strMaxCell = cell(lineNum,1);

yDataCompleteRange = [+inf,-inf];

for lineIdx = 1:lineNum
   yDataLineName{lineIdx} = strrep(strrep(tmpStrEval{lineIdx},'@',''),'_',' ');
   tmpStr = strrep(tmpStrEval{lineIdx},'@','ParametricResult.');
   yDataFullND = eval(tmpStr); % Compute full-ND matrix given z-expression

   % Find global max and min:
   [yDataGlobalMin, yDataGlobalMinIdx] = min(yDataFullND(:));
   [yDataGlobalMax, yDataGlobalMaxIdx] = max(yDataFullND(:));
   
   % Exchange global range if worse than previous line:
   if yDataCompleteRange(1) > yDataGlobalMin
      yDataCompleteRange(1) = yDataGlobalMin;
   end
   if yDataCompleteRange(2) < yDataGlobalMax
      yDataCompleteRange(2) = yDataGlobalMax;
   end

   tmpIdx = cell(1,length(size(yDataFullND)));
   [tmpIdx{:}] = ind2sub(size(yDataFullND),yDataGlobalMaxIdx);
   strMax = sprintf('max = %.3g at %s',yDataGlobalMax,...
      strrep(mat2str(cat(2,tmpIdx{:})),' ',','));
   tmpIdx = cell(1,length(size(yDataFullND)));
   [tmpIdx{:}] = ind2sub(size(yDataFullND),yDataGlobalMinIdx);
   strMin = sprintf('min = %.3g at %s',yDataGlobalMin,...
      strrep(mat2str(cat(2,tmpIdx{:})),' ',','));
   
   strMinCell{lineIdx} = strMin;
   strMaxCell{lineIdx} = strMax;
   yDataFullNDcell{lineIdx} = yDataFullND;   
end

% To make the selection arrays, the actual line does not matter:
activeDim = find(size(yDataFullND) > 1);
upToActiveDim = 1:max(activeDim);
permuteArray = unique(cat(2,xDim,upToActiveDim),'stable');
activeDim2control = setdiff(activeDim,[xDim]);

dimIsNonSingletonInExpr = setdiff(allDim,activeDim);
tmpSelectIdxCellArray = selectIdxCellArray;
tmpSelectIdxCellArray(dimIsNonSingletonInExpr) = {1};

% Create "plot" figure
hPlotFig = figure();

% Generate a control figure
hCtrFig           = figure('Name',sprintf('Plot control figure %d',hPlotFig.Number),...
   'Toolbar','none','Menubar','none','NumberTitle','off');

% Print uicointrol components start from the bottom of the list:
tmpCurrentPosition  = [GuiOpt.pixelPaddingWidth(1), GuiOpt.pixelPaddingHeigth(1)];

% Prepare small "control" window to the plotted figure:
for dimSelectIdx = length(activeDim2control):-1:1
   
   % Convert to actual dimension index:
   tmpDimIdx = activeDim2control(dimSelectIdx);
   
   % Select the vectorial data arrays we should display:
   tmpParameterPrint = parameterActiveDimIfVector == tmpDimIdx;
   
   % Create data arrays
   tmpParamStr  = parameterStrIfVector(tmpParameterPrint);
   tmpParamData = cellfun(@(str) ParametricResult.(str),tmpParamStr,'un',0);
   
   % Create "value" window
   tmpPosition   = [tmpCurrentPosition, ...
      GuiOpt.extentControlWindowTotal, ...
      length(tmpParamStr)*GuiOpt.pixelHeigthUIcontrol];
   
   tmpParamStrMaxLen = max(cellfun(@(str) length(str), tmpParamStr));
   tmpStrFormat       = ['%-' sprintf('%d',tmpParamStrMaxLen) 's = %s'];
   
   tmpStr  = cellfun(@(param,val) sprintf(tmpStrFormat, param, ...
      num2sci(val(selectedIdxInEachDim(tmpDimIdx)), 3)),...
      tmpParamStr,tmpParamData, 'un', 0);
   tmpStr = sprintf(cellstr2str(tmpStr,'\n'));
   
   hParText = uicontrol( ...
      'Style',       'edit',...
      'String',      tmpStr,...
      'Min',         0,...
      'Max',         2,...
      'Enable',      'inactive',...
      'FontName',    'fixedwidth', ...
      'FontSize',    GuiOpt.editFontSize, ...
      'Position',    tmpPosition, ...
      'HorizontalAlignment', 'left',...
      'Userdata',    {tmpParamStr, tmpParamData, tmpStrFormat});
   
   % Update current position
   tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
      length(tmpParamStr)*GuiOpt.pixelHeigthUIcontrol + ...
      GuiOpt.pixelPaddingHeigth(2);
   
   % Create "index text" field for slider
   tmpPosition   = [tmpCurrentPosition, ...
      GuiOpt.extentControlWindowText, ...
      GuiOpt.pixelHeigthUIcontrol];
   
   tmpStrFormat = 'dim %d, idx %d';
   tmpStr = sprintf(tmpStrFormat,tmpDimIdx,selectedIdxInEachDim(tmpDimIdx));
   hDimText = uicontrol( ...
      'Style',       'text',...
      'String',      tmpStr,...
      'FontName',    'times', ...
      'FontSize',    GuiOpt.textFontSize, ...
      'Position',    tmpPosition, ...
      'HorizontalAlignment', 'left',...
      'Userdata',    {tmpDimIdx, tmpStrFormat});
   
   % Update current position
   tmpCurrentPosition(1) = tmpCurrentPosition(1) + ...
      GuiOpt.extentControlWindowText + ...
      GuiOpt.pixelPaddingWidth(2);
   
   % Create "slider"
   tmpPosition   = [tmpCurrentPosition, ...
      GuiOpt.extentControlWindowSlider, ...
      GuiOpt.pixelHeigthUIcontrol];
   
   sliderStepsFraction = [1,5]./(parameterSizeAll(tmpDimIdx)-1);
   hSlider = uicontrol( ...
      'Style',       'slider',...
      'Min',         1,...
      'Max',         parameterSizeAll(tmpDimIdx),...
      'Value',       selectedIdxInEachDim(tmpDimIdx),...
      'Position',    tmpPosition, ...
      'Callback',    {@update_plot_slider, hDimText, hParText, hPlotFig},...
      'SliderStep',  sliderStepsFraction);
   
   % For continous updates when we move the slider we add a listener:
   if exist('addlistener','builtin')
      addlistener(hSlider,'Value','PostSet',@(s,e) update_plot_slider_listner(s,e,hDimText,hParText,hPlotFig));
   end
   
   % Update current position
   tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
      GuiOpt.pixelHeigthUIcontrol + ...
      GuiOpt.pixelPaddingHeigth(2);
   tmpCurrentPosition(1) = GuiOpt.pixelPaddingWidth(1);
   
end

% Print info window with min/max value and index:
tmpSize = 2.5*lineNum;
tmpPosition   = [tmpCurrentPosition, ...
   GuiOpt.extentControlWindowTotal, ...
   tmpSize*GuiOpt.pixelHeigthUIcontrol];

totalMinMaxStr = cellfun(@(strMin,strMax,lineIdx) ...
   sprintf('Line %d: maximum and minimum value\n%s\n%s\n',...
   lineIdx,strMin,strMax),strMinCell(:).',strMaxCell(:).',num2cell(1:lineNum),'un',0);
totalMinMaxStr = cat(2,totalMinMaxStr{:});

uicontrol( ...
   'Style',       'edit',...
   'String',      totalMinMaxStr,...
   'Min',         0,...
   'Max',         2,...
   'Enable',      'inactive',...
   'FontName',    'times', ...
   'FontSize',    GuiOpt.editFontSize, ...
   'Position',    tmpPosition, ...
   'HorizontalAlignment', 'left');

% Update current position
tmpCurrentPosition(2) = tmpCurrentPosition(2) + ...
   tmpSize*GuiOpt.pixelHeigthUIcontrol + ...
   GuiOpt.pixelPaddingHeigth(2);
   
% Update size of window
screenSize           = get(0,'ScreenSize');
figureSize           = [GuiOpt.extentControlWindowTotal + 2*GuiOpt.pixelPaddingWidth(1),...
   tmpCurrentPosition(2) + GuiOpt.pixelPaddingHeigth(1) - GuiOpt.pixelPaddingHeigth(2)];
figurePosition       = 0.5*screenSize([3,4]) - 0.5*figureSize;
set(hCtrFig,'Position',[figurePosition,figureSize]);

% Change -inf,+inf in zMaxMinRange to the actual max or min of the data set:
PlotOption.yMinMaxRange(PlotOption.yMinMaxRange == -inf) = yDataCompleteRange(1);
PlotOption.yMinMaxRange(PlotOption.yMinMaxRange == +inf) = yDataCompleteRange(2);

% Append important data to the userdata of the plot figure:
PlotUserData = struct();
PlotUserData.PlotOption = PlotOption;

PlotUserData.xData = xData;
% PlotUserData.yData = yData;
PlotUserData.yDataFullNDcell = yDataFullNDcell;
PlotUserData.selectedIdxInEachDim = selectedIdxInEachDim; % Currently selected idx
PlotUserData.yDataPermuteArray  = permuteArray;
PlotUserData.selectIdxCellArray = tmpSelectIdxCellArray;
PlotUserData.yDataLineName = yDataLineName;
PlotUserData.yDataCompleteRange = yDataCompleteRange;

set(hPlotFig,'Userdata',PlotUserData);

set(hPlotFig,'CloseRequestFcn',{@close_both_figures,hCtrFig});
set(hCtrFig,'CloseRequestFcn',{@close_both_figures,hPlotFig});
% Plot initial figure
plot_contour(hPlotFig,true)

% Update plot according to slider change
function update_plot_slider(hSlider,event,hDimText,hParText,hPlotFig)

% Get plot userdata structure:
PlotUserData = get(hPlotFig,'Userdata');
% PlotUserData.selectedIdxInEachDim;
% PlotUserData.selectIdxCellArray;

% Get dim of slider:
tmpUserData = get(hDimText,'Userdata');
tmpDim = tmpUserData{1};
tmpStrFormat = tmpUserData{2};

% Get value of slider:
tmpValue = get(hSlider,'Value');
tmpIntValue = round(tmpValue);
PlotUserData.selectedIdxInEachDim(tmpDim) = tmpIntValue;
PlotUserData.selectIdxCellArray{tmpDim} = tmpIntValue;

% Update text to the left of slider:
tmpStr = sprintf(tmpStrFormat, tmpDim, tmpIntValue);
set(hDimText,'String',tmpStr);

% Update values of parameters:
tmpUserData = get(hParText,'Userdata');
tmpParamStr = tmpUserData{1};
tmpParamData = tmpUserData{2};
tmpStrFormat = tmpUserData{3};
tmpStr  = cellfun(@(param,val) sprintf(tmpStrFormat, param, ...
   num2sci(val(PlotUserData.selectedIdxInEachDim(tmpDim)), 3)),...
   tmpParamStr,tmpParamData, 'un', 0);
tmpStr = sprintf(cellstr2str(tmpStr,'\n'));
set(hParText,'String',tmpStr);

% Update hPlotFig userdata
set(hPlotFig,'Userdata',PlotUserData)

% Update Plot:
plot_contour(hPlotFig)

% Update plot according to slider change
function update_plot_slider_listner(hListner,event,hDimText,hParText,hPlotFig)

% Fix for 2014b:
hSlider  = event.AffectedObject;

% Get plot userdata structure:
PlotUserData = get(hPlotFig,'Userdata');
% PlotUserData.selectedIdxInEachDim;
% PlotUserData.selectIdxCellArray;

% Get dim of slider:
tmpUserData = get(hDimText,'Userdata');
tmpDim = tmpUserData{1};
tmpStrFormat = tmpUserData{2};

% Get value of slider:
tmpValue = get(hSlider,'Value');
tmpIntValue = round(tmpValue);
PlotUserData.selectedIdxInEachDim(tmpDim) = tmpIntValue;
PlotUserData.selectIdxCellArray{tmpDim} = tmpIntValue;

% Update text to the left of slider:
tmpStr = sprintf(tmpStrFormat, tmpDim, tmpIntValue);
set(hDimText,'String',tmpStr);

% Update values of parameters:
tmpUserData = get(hParText,'Userdata');
tmpParamStr = tmpUserData{1};
tmpParamData = tmpUserData{2};
tmpStrFormat = tmpUserData{3};
tmpStr  = cellfun(@(param,val) sprintf(tmpStrFormat, param, ...
   num2sci(val(PlotUserData.selectedIdxInEachDim(tmpDim)), 3)),...
   tmpParamStr,tmpParamData, 'un', 0);
tmpStr = sprintf(cellstr2str(tmpStr,'\n'));
set(hParText,'String',tmpStr);

% Update hPlotFig userdata
set(hPlotFig,'Userdata',PlotUserData)

% Update Plot:
plot_contour(hPlotFig)
drawnow

% Function to plot contour
function plot_contour(hPlotFig,isInit)
if nargin==1
   isInit = false;
end

% Select the figure:
figure(hPlotFig)
if isInit
   box on
end
% Collect important data:
tmpUserData = get(hPlotFig,'Userdata');
PlotOption = tmpUserData.PlotOption;
xData = tmpUserData.xData;
% yData = tmpUserData.yData;
permuteArray = tmpUserData.yDataPermuteArray;

% if PlotOption.cfill
%    optFill = 'on';
% else
%    optFill = 'off';
% end

% Scale data according to dataScale:
xData = xData.*PlotOption.xScale;
% yData = yData.*PlotOption.yScale;
lineNum = length(tmpUserData.yDataFullNDcell);
for lineIdx = 1:lineNum
   yData = permute(tmpUserData.yDataFullNDcell{lineIdx}(...
      tmpUserData.selectIdxCellArray{:}),permuteArray);

% yDataRange  = [min(yData(:)), max(yData(:))];
% yDataCompleteRange = [min(tmpUserData.yDataFullND(:)),...
%    max(tmpUserData.yDataFullND(:))];

% if PlotOption.nLvlTotalRange
%    if ~isempty(PlotOption.yMinMaxRange)
%       useDataRange = PlotOption.yMinMaxRange;
%    else
%       useDataRange = yDataCompleteRange;
%    end
% else
%    useDataRange = yDataRange;
% end

% % Check for log-scale color:
% if PlotOption.zLogScale
%    funDoLog = @(x) log10(x);
% else
%    funDoLog = @(x) x;
% end

% Plot using contour:
% cLvl = PlotOption.nLvl;
% if isempty(cLvl)
%    cLvl = 6;
% end
% if isscalar(cLvl) && PlotOption.zLogScale
%    cLvl = logspace(log10(useDataRange(1)),log10(useDataRange(2)),cLvl);
% elseif isscalar(cLvl)
%    cLvl = linspace(useDataRange(1),useDataRange(2),cLvl);
% end
% 
% % Round levels
% cLvlFix = round(funDoLog(cLvl),3,'significant');
% 
% % Plot contour using contour or contourf:
% if PlotOption.cfill
%    [C,h] = contourf(xData,yData,funDoLog(yData),...
%       cLvlFix,'LineWidth',PlotOption.lw,...
%       'LineStyle',PlotOption.lt);
% else
%    [C,h] = contour(xData,yData,funDoLog(yData),...
%       cLvlFix,'LineWidth',PlotOption.lw,...
%       'LineStyle',PlotOption.lt);
% end

% Plot using plot:
plot(xData,yData,'LineWidth',PlotOption.lw,'LineStyle',PlotOption.lt)
if lineIdx == 1
   hold on
end
end
hold off



% Set x/y axis to logarithmic if specified:
if PlotOption.xLogScale
   set(gca,'XScale','log');
end
if PlotOption.yLogScale
   set(gca,'YScale','log');
end


% Set fontsize on axes:
grid on
set(gca,'fontsize',14)

% Set maximum y range:
if ~isempty(PlotOption.yMinMaxRange)
   xRange = [min(xData(:)), max(xData(:))];
   axis([xRange,PlotOption.yMinMaxRange])
else
   axis tight
end

% Print x,y labels and title
xlabel(PlotOption.xLabel)
ylabel(PlotOption.yLabel)
if ~isempty(PlotOption.title)
   title(PlotOption.title)
end

% Add legend if multiple lines
if lineNum > 1 && PlotOption.legend
   legend(tmpUserData.yDataLineName{:})
end

drawnow

function close_both_figures(object,event,hFig2)
delete(hFig2)
delete(object)