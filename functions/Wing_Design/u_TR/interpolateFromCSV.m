function value = interpolateFromCSV(files, xq, yq)
%valuta u, v, w 
%INput:
%pattern = l'inizio del nome del file csv: può essere AR*(u),TR*(v),ctcr*(w)
%xq,yq = punto query in cui voglio u,v,w 
files =dir('*.csv');
X = []; Y = []; Z = [];
    for k = 1:length(files)
        data = load(files(k).name);
        X = [X; data(:,1)];
        Y = [Y; data(:,2)];
        Z = [Z; data(:,3)];
    end
    F = scatteredInterpolant(X,Z,Y);
    value = F(xq,yq);
end