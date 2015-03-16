%% Parameter study generation script for elliptic integral 2 coil

%% Set up models to load.
inputModel = {'elliptic_two_coils_model'};
modelNum = length(inputModel);

%% Set up override parameters for input model
% Obs! @varName is a link to an existing parameters value.
% Note, if the same parameter name is used, the parameter is automatically the
% same. Therefore, it is not possible link Rp='@Rp' because it can not be
% overriden by the next model. If the variable name is new, it can link to any
% previously used input or output parameter via any vectorized expression.
% Obs! The @ sign is simply replaced with 'ParametricResults.', the expression
% is then evaluated using eval.

overrideParameterCell = cell(1, modelNum);

overrideParameterCell{1}.c1wire = 'l';
overrideParameterCell{1}.c2wire = 'l';
overrideParameterCell{1}.c1z0 = '-(@c1nz-1).*1e-2';
overrideParameterCell{1}.c2nr = '@c1nr';
overrideParameterCell{1}.c2nz = '@c1nz';
overrideParameterCell{1}.c2z0 = 3e-1;
overrideParameterCell{1}.c1r0 = '3e-1 - (@c1nr-1).*1e-2';
overrideParameterCell{1}.c2r0 = '3e-1 - (@c1nr-1).*1e-2';




%% Define parametric study. 
% Parameters subject to parametric sweep should be
% defined here. Example:
% exampleStruct.param1 = {dim, valueVector};
% Study sweeps all parameters with different dim. Parameters with equal
% dim updates together.
ParametricStudy = struct();

ParametricStudy.f = {1, linspace(70e3, 100e3, 50)};
ParametricStudy.c1nr = {2, 1:30};
ParametricStudy.c1nz = {3, 1:25};





%% Override parameterValues and generate ModelStruct input for parametric study:
modelStructCell = override_input_model_value(inputModel, overrideParameterCell);


%% Compute parametric study:
tic();
[ParametricResult, ParametricDimHelper] = compute_parametric_study(modelStructCell, ParametricStudy);
toc();

%% Store results
mkdirN('..\WETRICdata\analyt_2coil\')
[fileName,pathName,filterIndex] = uiputfile('.mat','Save parametric results','..\WETRICdata\analyt_2coil\');

save([pathName, fileName],'ParametricResult','modelStructCell')