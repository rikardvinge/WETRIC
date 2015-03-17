%% Parameter study generation script for analytic 2 coil, circuit and fields

%% Initially, we use two models
inputModel  = {'analyt_collapsed_2coil', 'circuit_2coil_parallel', 'field_amplitude_2coil'};
modelNum    = length(inputModel);

% Set up override values for the input models:
overrideParameter    = cell(1,modelNum);
% Override for input model 1:
overrideParameter{1}.N1 = 2;
overrideParameter{1}.N2 = 1;
overrideParameter{1}.d  = 0.30;
overrideParameter{1}.a1  = 0.30;
overrideParameter{1}.a2  = 0.30;
overrideParameter{1}.b  = 3e-3; 
overrideParameter{1}.Z  = 0.05; % Set z-coordinate
overrideParameter{1}.Y  = 0.00; % Set y-coordinate
% Override for input model 2:
overrideParameter{2}.ug = 1; % amplitude of voltage source
overrideParameter{2}.Rg = 0.25; % amplitude of voltage source
overrideParameter{2}.C1 = '1./((2*pi*85e3)^2.*@L1)';
overrideParameter{2}.C2 = '1./((2*pi*85e3)^2.*@L2)';
% Obs! @varName is a link to an existing parameters value.
% Note, if the same parameter name is used, the parameter is automatically the
% same. Therefore, it is not possible link R1='@R1' because it can not be
% overriden by the next model. If the variable name is new, it can link to any
% previously used input or output parameter via any vectorized expression.
% Obs! The @ sign is simply replaced with 'ParametricResults.', the expression
% is then evaluated using eval. Thus, if we use multiple parameters, it might be
% necessary to write the linked expressions using bsxfun.

%% Set parameters to do parametric study versus:
ParametricStudy      = struct();
% ParametricStudy.b    = {2, linspace(5e-4,5e-2,25)};
% ParametricStudy.N2   = {1, 1:1:20};
% ParametricStudy.d    = {2, linspace(0.10,0.30,11)};
ParametricStudy.f    = {1, linspace(80e3,90e3,101)};
ParametricStudy.Rl   = {2, linspace(5,100,20)};
ParametricStudy.X    = {3, 5*logspace(-2,1,25)};
% ParametricStudy.Y    = {5, 5*logspace(-2,1,20)};
% ParametricStudy.X    = {4, 5*logspace(-2,0,50)};
% ParametricStudy.Y    = {4, 5*logspace(-2,0,50)};
% ParametricStudy.Z    = {4, 5*logspace(-2,0,50)};
% ParametricStudy.Lp   = {4, linspace(10,60,51).*1e-6};
% ParametricStudy.Ls   = {4, linspace(10,60,51).*1e-6};
% Note, in a GUI this structure would be generated according to user input.

%% Set save name and folder:
mkdirN('..\WETRICdata\analyt_2coil\')
[fileName,pathName,filterIndex] = uiputfile('.mat','Save parametric results','..\WETRICdata\analyt_2coil\');

%% Override parameterValues and generate ModelStruct input for parametric study:
modelStructCell = override_input_model_value(inputModel, overrideParameter);

%% Compute parametric study:
[ParametricResult, ParametricDimHelper] = compute_parametric_study(modelStructCell, ParametricStudy);
% Note, I might remove the ParametricDimHelper output.

%% Store results
save([pathName, fileName],'ParametricResult','modelStructCell')