classdef SimpleApp < matlab.apps.AppBase

    % Componenti dell'app
    properties (Access = public)
        UIFigure    matlab.ui.Figure
        SalutaButton matlab.ui.control.Button
        ChiudiButton matlab.ui.control.Button
    end

    % Callbacks dei componenti
    methods (Access = private)

        % Callback del bottone "Saluta"
        function SalutaButtonPushed(app, event)
            disp('Ciao! Hai premuto il pulsante Saluta.');
        end

        % Callback del bottone "Chiudi"
        function ChiudiButtonPushed(app, event)
            disp('Chiusura dell''app.');
            delete(app.UIFigure);
        end
   

    % Codice di inizializzazione dei componenti

        % Crea e configura i componenti
        function createComponents(app)

            % Crea la finestra principale
            app.UIFigure = uifigure('Position', [500 500 300 150], ...
                                    'Name', 'Esempio GUI con App Designer');

            % Crea il bottone "Saluta"
            app.SalutaButton = uibutton(app.UIFigure, 'push');
            app.SalutaButton.Position = [50 80 200 40];
            app.SalutaButton.Text = 'Saluta';
            app.SalutaButton.ButtonPushedFcn = createCallbackFcn(app, @SalutaButtonPushed, true);

            % Crea il bottone "Chiudi"
            app.ChiudiButton = uibutton(app.UIFigure, 'push');
            app.ChiudiButton.Position = [50 20 200 40];
            app.ChiudiButton.Text = 'Chiudi';
            app.ChiudiButton.ButtonPushedFcn = createCallbackFcn(app, @ChiudiButtonPushed, true);

        end
    end

    % Costruttore
    methods (Access = public)

        % Costruttore dell'app
        function app = SimpleApp
            createComponents(app)
        end
    end
end