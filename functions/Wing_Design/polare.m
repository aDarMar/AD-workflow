function polare(CL,CD)
% Plotta un grafico in base ai vettori CL CD
% Input:
% CL = CL dell'ala 3d
% CD = CD dell'ala 3d

    % Verifica che i vettori abbiano la stessa lunghezza
    if length(CL) ~= length(CD)
        error('I vettori devono avere la stessa lunghezza.');
    end

    % Crea il grafico
    figure;
    polare(CL, CD, '-o', 'LineWidth', 2);
    xlabel('CL');
    ylabel('CD');
    title('CL/CD');
    grid on;

end