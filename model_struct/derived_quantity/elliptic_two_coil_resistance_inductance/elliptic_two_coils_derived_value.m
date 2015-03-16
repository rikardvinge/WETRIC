function [exprStr, plotLabel] = elliptic_two_coils_derived_value(expression)
% Generates expressiosn for derived values using bsxfun for singleton 
% extensions.
% Note: derived_value('implemented') returns information on implemented
% expressions.
%
% Main expressions:
% Simple expressions:
% k = M/sqrt(L1*L2)
% Q1 = wL1/R1
% Q2 = wL2/R2
% 
% Complex expression:
% QL_optim = Q2/sqrt(1 + k^2Q1Q2) -> [Kiani et al.]
% eta_2coil = k^2Q1Q2L/(1+k^2Q1Q2L) * Q2L/QL -> [Kiani et al.]
% Q1/w
% RL_optim = QL_optim*w*L2 -> [Kiani et al.]
%
% Secondary expressions (used in generation of main expressions):
% k^2*Q1
% k^2*Q1*Q2
% Q2L = Q2QL/(Q2+QL)
% k^2*Q1*Q2L
%
% Comment:
%
% Example usage:
% >> expressionString = derived_value('k');
% expressionString = 
% brdivide(@M, sqrt(btimes(@L1, @L2)))
%
% >> expressionString = derived_value('k');
% expressionString = 
% 
% See also bplus, bminus, btimes, brdivide

%   Created by: Rikard Vinge
%   $Revision: 1.0$  $Date: 2015-03-12 18:00:00$
%   $Revision: 1.01$ $Date: 2015-03-12 14:00:00$

% Control input:

if strcmp(expression, 'k')

  exprStr = 'brdivide(@M, sqrt(btimes(@L1,@L2)))';
  plotLabel = 'k [-]';
  
elseif strcmp(expression, 'Q1')
  exprStr = 'brdivide(btimes(2*pi*@f, @L1),@R1)';
  plotLabel = 'Q1 [-]';
  
elseif strcmp(expression, 'Q2')
  exprStr = 'brdivide(btimes(2*pi*@f, @L2),@R2)';
  plotLabel = 'Q2 [-]';
  
elseif strcmp(expression, 'QL_optim')
  % Optimal quality load factor according to Kiani et al.
  % QL_optim = Q2/sqrt(1 + k^2Q1Q2)
  k2Q1Q2Str = elliptic_two_coils_derived_value('k2Q1Q2');
  Q2Str = elliptic_two_coils_derived_value('Q2');
  denominatorStr = strcat('sqrt(1+',k2Q1Q2Str,')');
  QLStr = strcat('brdivide(',Q2Str,',',denominatorStr,')');
  exprStr = QLStr;
  plotLabel = 'QL [-]';
  
elseif strcmp(expression, 'eta_2coil')
  % Efficiency for 2-coil system according to Kiani et al.
  % eta = k^2Q1Q2/(1+k^2Q1Q2) * Q2L/QL
  QLStr = elliptic_two_coils_derived_value('QL_optim');
  Q2LStr = elliptic_two_coils_derived_value('Q2L');
  k2Q1Q2LStr = elliptic_two_coils_derived_value('k2Q1Q2L');
  firstFractionStr = strcat('brdivide(',k2Q1Q2LStr,',','1+',k2Q1Q2LStr,')');
  secondFractionStr = strcat('brdivide(',Q2LStr,',',QLStr,')');
  etaStr = strcat('btimes(',firstFractionStr,',',secondFractionStr,')');
  exprStr = etaStr;
  plotLabel = '\eta [-]';
  
elseif strcmp(expression, 'Q1/w')
  Q1Str = elliptic_two_coils_derived_value('Q1');
  exprStr = strcat('brdivide(',Q1Str,',2*pi*@f)');
  plotLabel = 'Q1/w [s]';
  
elseif strcmp(expression,'RL_optim')
  % RL_optim = QL_optim*w*L2 -> [Kiani et al.]
  QLStr = elliptic_two_coils_derived_value('QL_optim');
  wL2Str = 'btimes(2*pi*@f,@L2)';
  exprStr = strcat('btimes(',QLStr,',',wL2Str,')');
  plotLabel = 'RL [\Omega]';
  
  
  
  
  
  
elseif strcmp(expression,'k2Q1')
  kStr = elliptic_two_coils_derived_value('k');
  Q1Str = elliptic_two_coils_derived_value('Q1');
  k2Str = strcat(kStr,'.*',kStr);
  k2Q1Str = strcat('btimes(',k2Str,',',Q1Str,')');
  exprStr = k2Q1Str;
  plotLabel = '';
  
elseif strcmp(expression,'k2Q1Q2')
  Q2Str = elliptic_two_coils_derived_value('Q2');
  k2Q1Str = elliptic_two_coils_derived_value('k2Q1');
  k2Q1Q2Str = strcat('btimes(',k2Q1Str,',',Q2Str,')');
  exprStr = k2Q1Q2Str;
  plotLabel = '';
  
elseif strcmp(expression,'k2Q1Q2L')
  kStr = elliptic_two_coils_derived_value('k');
  Q1Str = elliptic_two_coils_derived_value('Q1');
  Q2LStr = elliptic_two_coils_derived_value('Q2L');
  k2Str = strcat(kStr,'.*',kStr);
  k2Q1Str = strcat('btimes(',k2Str,',',Q1Str,')');
  k2Q1Q2LStr = strcat('btimes(',k2Q1Str,',',Q2LStr,')');
  exprStr = k2Q1Q2LStr;
  plotLabel = '';
  
elseif strcmp(expression,'Q2L')
  % Loaded quality factor for 2-coil system according to Kiani et al.
  Q2Str = elliptic_two_coils_derived_value('Q2');
  QLStr = elliptic_two_coils_derived_value('QL_optim');
  denominatorStr = strcat(Q2Str,'+',QLStr);
  numeratorStr = strcat('btimes(',Q2Str,',',QLStr,')');
  Q2LStr = strcat('brdivide(',numeratorStr,',',denominatorStr,')');
  exprStr = Q2LStr;
  plotLabel = '';
  
else
  error(['Expression for ',expression,' not implemented.']);
  
end
