function [dCl0,cbaroc,alphaD] = dCl02DHLFun(df,flap,posFLG,Midx)%cloBase,a,cfoc,flaptype,posflag )
%dCl02D: Calcola la variazione di Cl0 su un profilo 2D a causa dei flap
%   Riceve in input la proprietà 

% K1 = k1FlapFun(cfc,flaptype,1);
% K2 = k2FlapFun(df,flaptype,1)
% K3 = k3FlapFun(dfodfref,flaptype,1)
% dClmaxBase = 

% tetaF = 1/cos(2*cfoc-1);
% alphaD = 1 -(tetaF - sin(tetaF))/pi;
% 
% etaD = etadelta_fFun(df,flaptype,1);%etaDeltaFun(df,1);
% dCoCf =  etaDeltaFlapFun(df,flaptype); % DEVE ESSERE FINITA
% cbaroc = 1+dCoCf*cfoc;
% dCl0bar = etaD*df*alphaD*a;
% 
% dCl0 = dCl0bar*cbaroc + cloBase*(cbaroc-1);

switch posFLG
    case 'tip'
        cloBase = flap.tip.cl0(Midx);
        a = flap.tip.a(Midx);
        cfoc = flap.tip.cfoc ;                                   %o,obj.flaps(i).
        flaptype = flap.tip.flaptype;
    case 'root'
        cloBase = flap.root.cl0(Midx);
        a = flap.root.a(Midx);
        cfoc = flap.root.cfoc ;                                   %o,obj.flaps(i).
        flaptype = flap.root.flaptype;
end

tetaF = 1/cos(2*cfoc-1);
alphaD = 1 -(tetaF - sin(tetaF/57.3))/pi;

etaD = etadelta_fFun(df,flaptype,1);%etaDeltaFun(df,1);
dCoCf =  etaDeltaFlapFun(df,flaptype); % DEVE ESSERE FINITA
cbaroc = 1+dCoCf*cfoc;

% switch posFLG 
%     case 'tip'
%     flap.tip.HLauxVariables('cbaroc',cbaroc);
%     case 'root'
%     flap.root.HLauxVariables('cbaroc',cbaroc);
% end

dCl0bar = etaD*df*alphaD*a;

dCl0 = dCl0bar*cbaroc + cloBase*(cbaroc-1);

end