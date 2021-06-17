clear
clc
close all

dx = 1; % cm
dy = 1; % cm
lambda = 1.75;% lambda de sobrerrelaxação
tolerancia = 0.0001; % criterio de parada do erro
limite_iters = 10000; % numero maximo de iterações
mi0 = 4*pi*10^-7;
miferro = 2500*mi0;
miar = mi0;
mibobina = mi0;

rows = 20/dy; % numero de linhas da matriz
cols = 22/dx; % numero de colunas da matriz
A = zeros(rows, cols); % matriz inicial com zeros
Anova = A;  % matriz com valores novos

% Calcula A por liebmann até alcançar a tolerância desejada ou bater o
% limite de iterações
while erro_max > tolerancia && iters < limite_iters
    
    iters = iters + 1;
    erro_max = 0;
        
    % Itera matriz por linha e coluna
    for j=1:rows-1 
        for i=1:cols-1
            
            pos_x = col*dx; % Posição do ponto no eixo x
            pos_y = row*dy; % Posição do ponto no eixo y
            
            % Se está embaixo ou em cima da bobina externa, pula
            if pos_x >= 20 && (pos_y <= 4 || pos_y >= 16)
                continue;
            end
                        
            % Fronteira vertical bobina-ferro esquerda
            if pos_x == 16 && pos_y > 4 && pos_y < 16
               
            % Fronteira vertical bobina-ferro direita
            elseif pos_x == 20 && pos_y > 4 && pos_y < 16        

            % Fronteira vertical ar-bobina
            elseif pos_x == 14 && pos_y > 4 && pos_y < 16 
            
            % Fronteira vertical ar-ferro esquerda
            elseif pos_x == 4
            
            % Fronteira vertical ar-ferro direita
            elseif (pos_y < 4 || pos_y > 16) && pos_x == 5
            
            % Fronteira horizontal bobina-ferro cima
            elseif pos_y == 16 && pos_x > 14 && pos_x < 16
            
            % Fronteira horizontal bobina-ferro baixo
            elseif pos_y == 4 && pos_x > 14 && pos_x < 16

            % Fronteira horizontal ar-ferro cima
            elseif pos_y == 16 && pos_x > 5 && pos_x < 14
            
            % Fronteira horizontal ar-ferro baixo
            elseif pos_y == 4 && pos_x > 5 && pos_x < 14
            
            % Quinas bobina-ferro-ar ???
            
            % Interior do domínio na bobina
            elseif pos_y > 4 && pos_y < 16 && ((pos_x > 14 && pos_x < 16) || pos_x > 20)      
                Acalc = (A(j, i+1) + A(j, i-1) + A(j+1, i) + A(j - 1, i))/4;
                
            % Interior do domínio ferro ou ar                      
            else
                Acalc = (A(j, i+1) + A(j, i-1) + A(j+1, i) + A(j - 1, i))/4 ...
                    + dx^2/4*(mi0*2*10^6*cos(pi*pos_y/(12*10^-2) + 8*10^5));
                
            end
            
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