clear
clc
close all

%% Parâmetros iniciais
[dx, dy] = deal(1); % passos em cm
lambda = 1.75;% lambda de sobrerrelaxação
tolerancia = 0.0001; % criterio de parada do erro
limite_iters = 10000; % numero maximo de iterações

rows = 20/dy; % numero de linhas da matriz
cols = 22/dx; % numero de colunas da matriz
A = zeros(rows, cols); % matriz inicial com zeros
Anova = A;  % matriz para guardar valores novos em uma iteração

mi0 = 4*pi*10^-7; % mi do vacuo
miferro = 2500*mi0; % mi do ferro
miar = mi0; % mi do ar
mibobina = mi0; % mi da bobina

%% Funções auxiliares para cálculo de Aij

% Função de cálculo de mi*Jz na bobina
miJz = @(y) mibobina*(2*10^6*cos(pi*y/(12*10^-2)) + 8*10^5); 

% Função de cálculo de Aij em um ponto no interior do domínio
Aij_interior = @(i, j, miJz)...
    (A(j, i+1) + A(j, i-1) + A(j+1, i) + A(j - 1, i))/4 + miJz;

% Função de cálculo de Aij na fronteira vertical entre dois meios
Aij_front_vert = @(i, j, mi1, Jz1, mi2, Jz2)...
    (-1/(dx*mibobina)*A(j, i-1) +  -1/(dx*mibobina)*A(j, i+1) ); 

% Função de cálculo de Aij na fronteira horizontal entre dois meios
Aij_front_hori = @(i, j, mi1, Jz1, mi2, Jz2) 1;
 
%% Aplicação do MDF

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
                
                Acalc = ((-1/(dx*mibobina)*A(j, i-1) +  (-1/(dx*mibobina)*A(j, i+1));
                
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
                Acalc = Aij_interior(i, j, miJz(pos_y));
                
            % Interior do domínio ferro ou ar                      
            else
                Acalc = Aij_interior(i, j, 0);
                
            end
            
            % Calcula valor novo de A no ponto por Liebmann
            Anova(j, i) = lambda*Acalc + (1-lambda)*A(j, i);

            % Calcula erro no ponto e substitui o erro máximo se for maior
            erro = abs(Anova(j, i) - A(j, i));
            if erro > erro_max
                erro_max = erro;
            end

        end
    end
    
    % Atualiza a matriz A
    A = Anova;
    
end