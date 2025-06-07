function [dClmax] = dClMaxHLFun(df,flap,posFLG,dref)%cfoc,toc,flaptype,dfodfref)
%dCl02D: Calcola la variazione di Cl0 su un profilo 2D a causa dei flap
%   Riceve in input la proprietà 

if nargin == 3
    dref = 40;
end

switch posFLG
    case 'tip'
        toc = flap.tip.tc;
        cfoc = flap.tip.cfoc ;                                   %o,obj.flaps(i).
        flaptype = flap.tip.flaptype;
    case 'root'
        toc = flap.root.tc;
        cfoc = flap.root.cfoc ;                                   %o,obj.flaps(i).
        flaptype = flap.root.flaptype;
end

K1 = k1FlapFun(cfoc,flaptype,1);
K2 = k2FlapFun(df,flaptype,1);
K3 = k3FlapFun(df/dref,flaptype,1);
dClmaxBase = dClMaxBaseFun(toc,flaptype,1);

dClmax = K1*K2*K3*dClmaxBase;

end