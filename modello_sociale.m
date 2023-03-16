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
prompt = "Si considerano interazioni asimmetriche o simmetriche? (si inserica 0 per asimmetriche e 1 per simmetriche) ";
simmetriche = input(prompt); % assume 1 per interazioni simmetriche, 0 per interazioni asimmetriche
prompt = "Che genere di distribuzione iniziale è presente? (1 - equispaziato; 2 - random; 3 - normale) ";
input_pos = input(prompt); % assume 1 per equispaziato, 2 per random, 3 per distribuzione normale
if input_pos == 3
    prompt = "Con che deviazione standard? (numero tra 0 e 1) ";
    deviazione = input(prompt); % deviazione standard
end

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
            f = f + x(i) * A(:,j) - x(j) * A(:,j);
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
            f = f + x(i) * A(:,j) - x(j) * A(:,j);
        end
        xtempo(:, indice) = x + dt * f; % applicazione di eulero
    end
end

for i=1:n
    figure(1)
    plot(t,xtempo(i,:));
    hold on
end