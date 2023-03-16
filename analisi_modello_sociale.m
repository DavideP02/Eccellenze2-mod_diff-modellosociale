% Questo file è identico al file "modello_sociale", ma si cerca di
% analizzare, al variare dei parametri iniziali, i risultati

clear all; close all; clc; % pulisce un po' tutto (variabili, finestre aperte, 
% command windows)

dt = 0.1; % fisso il passo sul tempo
Tmax = 20; % fisso la fine del processo
t = 0:dt:Tmax; % creo il vettore tempo
colonne = length(t); % numero di punti temporali 

% input
prompt = "Quanti agenti sono in gioco? ";
n = input(prompt);
prompt = "Quale distanza di influenza si considera? ";
epsilon = input(prompt);

% valori che variano:
simmetriche = 1; % assume 1 per interazioni simmetriche, 0 per interazioni asimmetriche
input_pos = 1; % assume 1 per equispaziato, 2 per random, 3 per distribuzione normale
deviazione = 0.2; % deviazione standard

titoli = ['Simmetrico - Equispaziato', 'Simmetrico - Random', 'Simmetrico - Normale', 'Asimmetrico - Equispaziato', 'Asimmetrico - Random', 'Asimmetrico - Normale'];


for situazioneiniziale = 1:6
    
    % l'indice situazioneiniziale è sul vettore titoli. Modifico i
    % parametri affinché coincidano
    if situazioneiniziale < 4
        simmetriche = 1;
    elseif situazioneiniziale > 3
        simmetriche = 0;
    end
    if situazioneiniziale == 1 || situazioneiniziale == 4 % equispaziato
        input_pos = 1;
    elseif situazioneiniziale == 2 || situazioneiniziale == 5 % random
        input_pos = 2; 
    elseif situazioneiniziale == 3 || situazioneiniziale == 6 % normale
        input_pos = 3;
    end

    if input_pos == 3
        for deviazione=0:0.05:1
            % vettore delle posizioni iniziali
            x0 = zeros(n,1); % definisco un vettore colonna di n uni
            if input_pos == 1
                % POSIZIONE INIZIALE EQUISPAZIATA
                h = 1 / (n - 1); % passo
                x0(1) = 0;
                for i=2:n
                    x0(i) = x0(i-1) + h;
                end
            elseif input_pos == 2
                % POSIZIONE INIZIALE RANDOM
                x0 = rand(n,1);
            elseif input_pos == 3
                % POSIZIONE INIZIALE DISTRIBUZIONE NORMALE
                x0 = 0.5 + deviazione .* randn(n,1);
            end
            
            % vettore delle posizioni
            xtempo = zeros(n,colonne); % inizializzazione
            xtempo(:,1) = x0; % gli assegno i valori iniziali
            
            % applicazione del metodo di EE
            
            for indice=2:colonne
                % adesso lavoro con indice-1 per calcolare f(t(indice-1), x(indice-1))
                x = xtempo(:, indice - 1); % creo il vettore posizione al tempo indice - 1
                if simmetriche == 1
                    % creo la matrice A
                    A = zeros(n,n);
                    for i=1:n
                        for j=1:n
                            if abs(x(i)-x(j)) < epsilon
                                A(i,j) = 1 / n;
                            end
                        end
                    end 
                    % applico il metodo
                    f = zeros(n,1); % calcolo la funzione per il tempo (indice - 1)
                    for j=1:n
                        f = f + x(j) * A(:,j) - x(i) * A(:,j);
                    end
                    xtempo(:, indice) = x + dt * f; % applicazione di eulero
                elseif simmetriche == 0
                    % devo calcolare ora i vari sigma(i)
                    sigma = zeros(n,1);
                    for i=1:n
                        for j=1:n
                            if abs(x(i)-x(j)) < epsilon
                                sigma(i) = sigma(i) + 1;
                            end
                        end
                    end
                    % creo la matrice A
                    A = zeros(n,n);
                    for i=1:n
                        for j=1:n
                            if abs(x(i)-x(j)) < epsilon
                                A(i,j) = 1 / sigma(i);
                            end
                        end
                    end
                    % applico il metodo
                    f = zeros(n,1); % calcolo la funzione per il tempo (indice - 1)
                    for j=1:n
                        f = f + x(j) * A(:,j) - x(i) * A(:,j);
                    end
                    xtempo(:, indice) = x + dt * f; % applicazione di eulero
                end
            end
            
            variazioni = zeros(1, colonne);
            for indice = 1:colonne
                variazioni(indice) = abs(max(xtempo(:,indice)) - min(xtempo(:,indice)));
            end
            
                    figure(situazioneiniziale)
        plot(t,variazioni);
        if situazioneiniziale == 1
            title('Simmetrico - Equispaziato')
        elseif situazioneiniziale == 2
            title('Simmetrico - Random')
        elseif situazioneiniziale == 3
            title('Simmetrico - Distribuzione Normale')
        elseif situazioneiniziale == 4
            title('Asimmetrico - Equispaziato')
        elseif situazioneiniziale == 5
            title('Asimmetrico - Random')
        elseif situazioneiniziale == 6
            title('Asimmetrico - Distribuzione Normale')
        end
        hold on
        end
    else
        % vettore delle posizioni iniziali
        x0 = zeros(n,1); % definisco un vettore colonna di n uni
        if input_pos == 1
            % POSIZIONE INIZIALE EQUISPAZIATA
            h = 1 / (n - 1); % passo
            x0(1) = 0;
            for i=2:n
                x0(i) = x0(i-1) + h;
            end
        elseif input_pos == 2
            % POSIZIONE INIZIALE RANDOM
            x0 = rand(n,1);
        elseif input_pos == 3
            % POSIZIONE INIZIALE DISTRIBUZIONE NORMALE
            x0 = 0.5 + deviazione .* randn(n,1);
        end
        
        % vettore delle posizioni
        xtempo = zeros(n,colonne); % inizializzazione
        xtempo(:,1) = x0; % gli assegno i valori iniziali
        
        % applicazione del metodo di EE
        
        for indice=2:colonne
            % adesso lavoro con indice-1 per calcolare f(t(indice-1), x(indice-1))
            x = xtempo(:, indice - 1); % creo il vettore posizione al tempo indice - 1
            if simmetriche == 1
                % creo la matrice A
                A = zeros(n,n);
                for i=1:n
                    for j=1:n
                        if abs(x(i)-x(j)) < epsilon
                            A(i,j) = 1 / n;
                        end
                    end
                end 
                % applico il metodo
                f = zeros(n,1); % calcolo la funzione per il tempo (indice - 1)
                for j=1:n
                    f = f + x(j) * A(:,j) - x(i) * A(:,j);
                end
                xtempo(:, indice) = x + dt * f; % applicazione di eulero
            elseif simmetriche == 0
                % devo calcolare ora i vari sigma(i)
                sigma = zeros(n,1);
                for i=1:n
                    for j=1:n
                        if abs(x(i)-x(j)) < epsilon
                            sigma(i) = sigma(i) + 1;
                        end
                    end
                end
                % creo la matrice A
                A = zeros(n,n);
                for i=1:n
                    for j=1:n
                        if abs(x(i)-x(j)) < epsilon
                            A(i,j) = 1 / sigma(i);
                        end
                    end
                end
                % applico il metodo
                f = zeros(n,1); % calcolo la funzione per il tempo (indice - 1)
                for j=1:n
                    f = f + x(j) * A(:,j) - x(i) * A(:,j);
                end
                xtempo(:, indice) = x + dt * f; % applicazione di eulero
            end
        end
        
        variazioni = zeros(1, colonne);
        for indice = 1:colonne
            variazioni(indice) = abs(max(xtempo(:,indice)) - min(xtempo(:,indice)));
        end

        titoli = ['Simmetrico - Equispaziato', 'Simmetrico - Random', 'Simmetrico - Distribuzione Normale', 'Asimmetrico - Equispaziato', 'Asimmetrico - Random', 'Asimmetrico - Distribuzione Normale'];
        
        figure(situazioneiniziale)
        plot(t,variazioni);
        if situazioneiniziale == 1
            title('Simmetrico - Equispaziato')
        elseif situazioneiniziale == 2
            title('Simmetrico - Random')
        elseif situazioneiniziale == 3
            title('Simmetrico - Distribuzione Normale')
        elseif situazioneiniziale == 4
            title('Asimmetrico - Equispaziato')
        elseif situazioneiniziale == 5
            title('Asimmetrico - Random')
        elseif situazioneiniziale == 6
            title('Asimmetrico - Distribuzione Normale')
        end
        hold on
        
    end
end