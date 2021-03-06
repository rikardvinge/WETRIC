%% Symbolic derivation of modified circuit model krets #5, 2 coils, parallel C1 
% Circuit model
%  (1) krets #5 v2 (without "old" C2)
%               vC1                   v4
%      o--[Rg]--o--------o   o--------o------o
%      |   <iRg |  iL1>  |   |  <iL2  | iRl> |
%    (ug)      C1\/iC1  L1:M:L2      C2\/iC2 Rl
%      |        |        |   |        |      |
%      o--------o--[R1]--o   o--[R2]--o------o
%              GND      vR1  v5        v6

%% Set up descrition of input and output symbols and desc. of circuit design:
CircDesc                   = struct();
CircDesc.functionName      = {'circuit_2coil_parallel','circuit_2coil_parallel_derived_quantites'};
CircDesc.functionSupportArrayInput = [true,true];
CircDesc.inputParameter    = {{'ug'},...
   {'Rg','C1','L1','R1','M','L2','R2','C2','Rl','f'},...
   {'iRg','iL1','iC1','iL2','iC2','iRl','vRl','vRg','vC1','vR1','vR2'}};
f0 = 85e3;
w0 = f0.*(2*pi);
k  = 0.10;
L  = 30e-6;
M  = k*L;
CircDesc.inputValue        = {{1},...
   {0.25,1/(w0^2*L),L,3e-3,M,L,3e-3,1/(w0^2*L),5,f0},...
    1,1,1,1,1,1,1,1,1,1,1};
CircDesc.outputParameter   = {...
   {'iRg','iL1','iC1','iL2','iC2','iRl','vRl','vRg','vC1','vR1','vR2'},...
   {'Pug','Pin','Pl','P1','P2','Prg','eta','eta_R'}};
% CircDesc.outputValue       = {};

CircDesc.design  = {...
'         vC1                     v4';...
'  o--[Rg]--o--------o   o--------o------o';...
'  |   <iRg |  iL1>  |   |  <iL2  | iRl> |';...
'(ug)      C1\/iC1  L1:M:L2      C2\/iC2 Rl';...
'  |        |        |   |        |      |';...
'  o--------o--[R1]--o   o--[R2]--o------o';...
'          GND      vR1  v5       v6';...
};
CircDesc.designDesc = {...  
'vRg = ug - vC1';...
'vR2 = v5 - v6';...
'vRl = v4 - v6';...
'vC2 = vRl';...
};

%% Set up circuit model [krets #5 v2]:
syms Rg C1 L1 R1 M L2 R2 C2 Rl f positive
syms ug real
piS = sym(pi);
w = 2*piS*f;

% Solve linear system of equations:
aMtx = sym(zeros(4));
aMtx(1,1) = 1/Rg + 1i*w*C1;
aMtx(1,2) = 1;
aMtx(2,1) = -1;
aMtx(2,2) = R1 + 1i*w*L1; % + 1/(1i*w*C2);
aMtx(2,4) = 1i*w*M;
aMtx(3,3) = 1/Rl + 1i*w*C2;
aMtx(3,4) = 1;
aMtx(4,2) = 1i*w*M;
aMtx(4,3) = -1;
aMtx(4,4) = R2 + 1i*w*L2;    
bVtr = sym(zeros(4,1));
bVtr(1) = ug/Rg;
xVtr = simplify(aMtx\bVtr);
% Derived quantities:
vC1 = xVtr(1);
iL1 = xVtr(2);
vRl = xVtr(3);
iL2 = xVtr(4);

iRg = (vC1-ug)/Rg;
iC1 = 1i*w*C1*vC1;
iC2 = 1i*w*C2*vRl;
iRl = vRl/Rl;

% v2 = vC1; % - i3/(1i*2*pi*f*C2);
vR1 = R1*iL1;
vR2 = R2*iL2;
vRg = ug - vC1;
    
% Generate Matlab function from symbolic expressions:
tmpAllOutSym = {iRg,iL1,iC1,iL2,iC2,iRl,vRl,vRg,vC1,vR1,vR2};


matlabFunction(tmpAllOutSym{:},'vars',cat(2,CircDesc.inputParameter{[1,2]}),...
   'outputs',CircDesc.outputParameter{1},...
   'file',[pwd '\model_struct\circuit\' CircDesc.functionName{1}]);
save([pwd '\model_struct\circuit\' CircDesc.functionName{1}],'-struct','CircDesc')


% tmpAllOutFunc = cell(1,length(tmpAllOutSym));
% for outIdx = 1:length(tmpAllOutSym)
%    tmpAllOutFunc{outIdx} = matlabFunction(tmpAllOutSym{outIdx},...,
%       'vars',CircDesc.inputParameter,...
%       'outputs',CircDesc.outputParameter);
% end

%% Set up derived quantities:

% Express formulas in previous functions output parameters:
syms ug iRg iL1 iC1 iL2 iC2 iRl vRl vRg vC1 vR1 vR2

% Power:
Pug   = ug*conj(-iRg);     % generated by voltage source
Pin   = vC1*conj(-iRg);    % input power to resonant circuit [remove loss in Rg]
Pl    = vRl*conj(iRl);     % load power
% Resistors:
P1    = vR1*conj(iL1);
P2    = vR2*conj(iL2);
Prg   = vRg*conj(-iRg);
% Efficiency: % %'real(@vRl.*conj(@iRl)./(@vC1.*conj(-@iRg)))';
eta   = simplify(Pl/Pin);         % input power versus power to load
PtotR = P1 + P2 + Pl + Prg;
eta_R = simplify(Pl/PtotR);    % power in load resistance versus power in all resistors


tmpAllOutSym = {Pug,Pin,Pl,P1,P2,Prg,eta,eta_R};
matlabFunction(tmpAllOutSym{:},'vars',cat(2,CircDesc.inputParameter{[1,3]}),...
   'outputs',CircDesc.outputParameter{2},...
   'file',[pwd '\model_struct\circuit\' CircDesc.functionName{2}]);