function [R1, R2] = elliptic_two_coils_resistance(f, ...
  nr1, nz1, dr1, r01, a1, sigma1, wireType1, ...
  nr2, nz2, dr2, r02, a2, sigma2, wireType2)
% Calculates the resistance of nd-matrix of geometry for one single coil.
% Note that wireType must be single-valued cell, it cannot contain multiple
% values.

% Constants.
litzWireDensity = 0.9069; % 0.9069 = Optimal packing density of circles in circle.
mu0 = 4*pi*1e-7;

Rtmp = cell(2,1);

for iTmp = 1:2
  if (iTmp == 1)
    nr = nr1;
    nz = nz1;
    dr = dr1; 
    r0 = r01;
    a = a1;
    sigma = sigma1;
    wireType = wireType1;
  else
    nr = nr2;
    nz = nz2;
    dr = dr2; 
    r0 = r02;
    a = a2;
    sigma = sigma2;
    wireType = wireType2;
  end
  % Length of coil wire, calculated by property of the aritmetic sum.
  wireLength = 2.*pi.*nz.*nr.*(r0 + 0.5.*(nr-1).*dr);
  
  % Calculate effective area passing current.
  if strcmp(wireType, 's');
    deltaTmp = sqrt(2./(2.*pi.*mu0.*f.*sigma));
    deltaTmp(10.*deltaTmp(:) > a(:)) = nan;
    area = pi.*a.*deltaTmp;
  elseif strcmp(wireType, 't');
    area = pi.*a.*a;
  elseif strcmp(wireType, 'l');
    area = pi.*a.*a.*litzWireDensity;
  else
    error('Must specify coil wire type as ''skin'', ''thin'', or ''litz''.');
  end
  
  % Calculate resistance.
  Rtmp{iTmp} = wireLength./(sigma.*area);
end

R1 = Rtmp{1};
R2 = Rtmp{2};

end