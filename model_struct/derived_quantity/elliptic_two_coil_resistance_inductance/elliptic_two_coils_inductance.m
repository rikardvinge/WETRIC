function [L1, L2, M] = elliptic_two_coils_inductance(...
  c1nr, c1nz, c1r0, c1z0, c1dr, c1dz, c1a,...
  c2nr, c2nz, c2r0, c2z0, c2dr, c2dz, c2a)
% Calculates the self- and mutual inductance of two coils for one set of
% coils per call.
% Algorithm based on Thomas Rylander's code.

mu0 = 4.*pi.*1e-7;


nCoils = 2;
rCell = cell(nCoils,1);
zCell = cell(nCoils,1);
% nrCell = cell(nCoils,1);
% nzCell = cell(nCoils,1);
aCell = cell(nCoils,1);

[c1rMtx, c1zMtx] = ndgrid((c1r0+((1:c1nr)-1)*c1dr), (c1z0+((1:c1nz)-1)*c1dz));
[c2rMtx, c2zMtx] = ndgrid((c2r0+((1:c2nr)-1)*c2dr), (c2z0+((1:c2nz)-1)*c2dz));
rCell{1} = c1rMtx(:);
rCell{2} = c2rMtx(:);
zCell{1} = c1zMtx(:);
zCell{2} = c2zMtx(:);
drCell{1} = c1dr;
drCell{2} = c2dr;
dzCell{1} = c1dz;
dzCell{2} = c2dz;
% nrCell{1} = c1nr;
% nrCell{2} = c2nr;
% nzCell{1} = c1nz;
% nzCell{2} = c2nz;
aCell{1} = c1a;
aCell{2} = c2a;

LMtx = zeros(nCoils, nCoils);

for traIdx = 1:nCoils
  
  % Transmitter
  dr_tra = drCell{traIdx};
  dz_tra = dzCell{traIdx};
  r_tra = rCell{traIdx};
  z_tra = zCell{traIdx};
  a_tra = aCell{traIdx};
  
  if any((dr_tra(:) < 3*a_tra(:)) | (dz_tra(:) < 3*a_tra(:)))
    % Should this be changed to return NaN
    error('Proximity assumption violated')
  end
  
  
  for recIdx = traIdx:nCoils
    
    % --- reciever ---
    r_rec = rCell{recIdx}.';
    z_rec = zCell{recIdx}.';
    a_rec = aCell{recIdx};
    
    % Inductance by bsxfunc
    aTmp = r_rec;
    bTmp = r_tra;
    % Use bsxfun to calculate nd-arrays of mTmp.
    % Note order: row: bTmp, col:aTmp.
    % cTmp = |z_rec - z_tra|
    cTmp = (bsxfun(@minus, z_tra, z_rec));
    % a+b = r_rec + r_tra
    AplusB = bsxfun(@plus, aTmp, bTmp);
    % a*b = r_rec.*r_tra
    AtimesB = bsxfun(@times, aTmp, bTmp);
    % m = 4*(a*b)./((a+b).^2 + c.^2)
    mTmp = 4.*AtimesB./((AplusB).^2 + cTmp.^2);
    
    
    % Calculate elliptical integral
    [kEiTmp, eEiTmp] = ellipke(mTmp);
    
    dEiTmp = (kEiTmp-eEiTmp)./mTmp;
    bEiTmp = (kEiTmp-dEiTmp);
    cEiTmp = (dEiTmp-bEiTmp)./mTmp;
    
    % If receiver == transiver: Swap diagonal values to self-inductance.
    if (traIdx == recIdx)
      [nRow, nCol] = size(cEiTmp);
      diagonalIdx = 1:(nRow+1):nRow*nCol;
      cEiTmp(diagonalIdx) = log(8.*r_rec(:)./a_rec) - 2;
    else
      if max(mTmp(:)) > 1%0.9999
        keyboard
      end
    end
    resTmp = sum(mu0.*sqrt(AtimesB(:)).*(mTmp(:).^(3/2)).*cEiTmp(:));
    LMtx(traIdx, recIdx) = resTmp;
    LMtx(recIdx, traIdx) = resTmp;
    
  end
end

L1 = LMtx(1,1);
L2 = LMtx(2,2);
M = LMtx(1,2);
end