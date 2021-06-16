clear
clc
close all

dx = 1; % cm
dy = 1; % cm
lambda = 1.75;% lambda de sobrerrelaxação
tolerancia = 0.0001; % criterio de parada do erro
limite_iters = 10000; % numero maximo de iterações

rows = 20/dy; % numero de linhas da matriz
cols = 22/dx; % numero de colunas da matriz
A = zeros(rows, cols); % matriz inicial com zeros
Anova = A;  % matriz com valores novos

colfinal = cols-1; % ultima coluna que será calculada na matriz (contorno = 0)
linha_inicio_bobina = 6/20*rows; % linha em que começa a bobina
linha_fim_bobina = 18/20*rows; % linha em que termina a bobina

% Calcula A por liebmann até alcançar a tolerância desejada ou bater o
% limite de iterações
while erro_max > tolerancia && iters < limite_iters
    
    iters = iters + 1;
    erro_max = 0;

    % Itera matriz por linha e coluna
    for row=1:rows-1 
        for col=1:cols-1
            % Se está embaixo ou em cima da bobina externa, pula
            if col == colfinal && (row < linha_inicio_bobina || row > linha_fim_bobina)
                continue;
            end
            
            pos_x = row*dx;
            pos_y = row*dy;
            
            % Calcula valor do ponto
            Acalc = formuleta;

            % Interior do domínio ferro ou ar
            
            % Interior do domínio na bobina

            % Fronteira vertical bobina-ferro esquerda
            
            % Fronteira vertical bobina-ferro direita

            % Fronteira vertical ar-bobina

            % Fronteira vertical ar-ferro esquerda
            
            % Fronteira vertical ar-ferro direita

            % Fronteira horizontal bobina-ferro cima
            
            % Fronteira horizontal bobina-ferro baixo

            % Fronteira horizontal ar-ferro cima
            
            % Fronteira horizontal ar-ferro baixo

            % Calcula valor novo de A no ponto por Liebmann
            Anova(row, col) = lambda*Acalc + (1-lambda)*A(row,col);

            % Calcula erro no ponto e substitui o erro máximo se for maior
            erro = abs(Anova(row, col) - A(row, col));
            if erro > erro_max
                erro_max = erro;
            end

        end
    end

    A = Anova;
    
end