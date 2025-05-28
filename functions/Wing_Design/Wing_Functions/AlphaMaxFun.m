function [alphaMax] = AlphaMaxFun(clmax,a,alpha0l,dY,sweep)
%AlphaMaxFun Funzione che calcola l0 alpha di stallo dell'ala 3D
dAlphaMax = dAlphaMaxFun(sweep,dY);
alphaMax = clmax/a + alpha0l + dAlphaMax;
end

