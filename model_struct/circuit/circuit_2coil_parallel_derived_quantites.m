function [Pug,Pin,Pl,P1,P2,Prg,eta,eta_R] = circuit_2coil_parallel_derived_quantites(ug,iRg,iL1,iC1,iL2,iC2,iRl,vRl,vRg,vC1,vR1,vR2)
%CIRCUIT_2COIL_PARALLEL_DERIVED_QUANTITES
%    [PUG,PIN,PL,P1,P2,PRG,ETA,ETA_R] = CIRCUIT_2COIL_PARALLEL_DERIVED_QUANTITES(UG,IRG,IL1,IC1,IL2,IC2,IRL,VRL,VRG,VC1,VR1,VR2)

%    This function was generated by the Symbolic Math Toolbox version 6.1.
%    17-Mar-2015 10:17:20

t2 = conj(iRg);
Pug = -t2.*ug;
if nargout > 1
    Pin = -t2.*vC1;
end
if nargout > 2
    t3 = conj(iRl);
    t4 = conj(iL1);
    t5 = t4.*vR1;
    t6 = conj(iL2);
    t7 = t6.*vR2;
    t8 = t3.*vRl;
    Pl = t8;
end
if nargout > 3
    P1 = t5;
end
if nargout > 4
    P2 = t7;
end
if nargout > 5
    Prg = -t2.*vRg;
end
if nargout > 6
    eta = -(t3.*vRl)./(t2.*vC1);
end
if nargout > 7
    eta_R = (t3.*vRl)./(t5+t7+t8-t2.*vRg);
end
