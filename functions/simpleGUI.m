function simpleGUI()
    % Crea la figura
    f = figure('Position', [100, 100, 800, 400], 'Name', 'Interfaccia Grafica con Grafici');

    % Campo di input per il parametro 'a'
    uicontrol('Style', 'text', 'Position', [50, 350, 100, 20], ...
              'String', 'Valore a:', 'HorizontalAlignment', 'left');
    a_input = uicontrol('Style', 'edit', 'Position', [150, 350, 100, 25], ...
                        'String', '1');

    % Pulsante per aggiornare i grafici
    uicontrol('Style', 'pushbutton', 'Position', [270, 350, 100, 25], ...
              'String', 'Aggiorna', 'Callback', @update_plots);

    % Assi per i due grafici
    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);

    % Inizializza i grafici
    x = linspace(0, 2*pi, 1000);
    a = str2double(a_input.String);
    h1 = plot(ax1, x, a*sin(x));
    title(ax1, 'a*sin(x)');
    h2 = plot(ax2, x, a*cos(x));
    title(ax2, 'a*cos(x)');

    % Funzione di callback per aggiornare i grafici
    function update_plots(~, ~)
        a = str2double(a_input.String);
        if isnan(a)
            errordlg('Inserire un valore numerico valido per a.');
            return;
        end
        h1.YData = a * sin(x);
        h2.YData = a * cos(x);
    end
end

