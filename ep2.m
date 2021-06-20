clear
clc
close all

%% Parâmetros iniciais
[dx, dy] = deal(0.01); % passos em metros
lambda = 1.75;% lambda de sobrerrelaxação
tolerancia = 0.0001; % criterio de parada do erro
limite_iters = 10000; % numero maximo de iterações

rows = 0.21/dy; % numero de linhas da matriz
cols = 0.23/dx; % numero de colunas da matriz
A = zeros(rows, cols); % matriz inicial com zeros

% Anova = A;  % matriz para guardar valores novos em uma iteração

mi0 = 4*pi*10^-7; % mi do vacuo
miferro = 2500*mi0; % mi do ferro
miar = mi0; % mi do ar
mibobina = mi0; % mi da bobina

% mi = ones(rows, cols)*NaN;
% mi(:, 1:4/23*cols) = miferro;
% mi(1:4/23*cols, :) = miferro;

%% Funções auxiliares para cálculo de Aij

% Função de cálculo de Jz na bobina
Jz = @(y) (2*10^6*cos(pi*y/(12*10^-2)) + 8*10^5); 

% Função de cálculo de Aij em um ponto no interior do domínio
Aij_interior = @(A, j, i, mi, Jz)...
    (A(j, i+1) + A(j, i-1) + A(j+1, i) + A(j - 1, i))/4 + mi*Jz;

% Função de cálculo de Aij na fronteira vertical entre dois meios
Aij_front_vert = @(A, j, i, mi1, Jz1, mi2, Jz2)...
    ((-1/(dx*mi1))*A(j, i-1) + (-1/(dx*mi2))*A(j, i+1) + ...
    (-1/(2*dx))*(1/mi2 + 1/mi1)*(A(j+1, i) + A(j-1, i)) + ...
    (-dx/2)*(Jz2 + Jz1))/((1/mi2 + 1/mi1)*(-2/dx)); 

% Função de cálculo de Aij na fronteira horizontal entre dois meios
Aij_front_hori = @(A, j, i, mi1, Jz1, mi2, Jz2) ...
    ((-1/(dx*mi1))*A(j-1, i) + (-1/(dx*mi2))*A(j+1, i) + ...
    (-1/(2*dx))*(1/mi2 + 1/mi1)*(A(j, i+1) + A(j, i-1)) + ...
    (-dx/2)*(Jz2 + Jz1))/((1/mi2 + 1/mi1)*(-2/dx));
 
%% Aplicação do MDF

% Calcula A por liebmann até alcançar a tolerância desejada ou bater o
% limite de iterações
iters = 0;
erro_max = inf;
while erro_max > tolerancia && iters < limite_iters
    
    iters = iters + 1;
    erro_max = 0;
    Aold = A;
    
    % Itera matriz por linha e coluna
    for j=2:rows-1 
        for i=2:cols-1
            
            pos_x = (i-1)*dx; % Posição do ponto no eixo x
            pos_y = (j-1)*dy; % Posição do ponto no eixo y
            
            % Se está embaixo ou em cima da bobina externa, pula
            if pos_x >= 0.20 && (pos_y <= 0.04 || pos_y >= 0.16)
                continue;
            end
                        
            % Fronteira vertical bobina-ferro esquerda
            if pos_x == 0.16 && pos_y > 0.04 && pos_y < 0.16
                Acalc = Aij_front_vert(A, j, i, mibobina, Jz(pos_y), miferro, 0);
            
            % Fronteira vertical bobina-ferro direita
            elseif pos_x == 0.20 && pos_y > 0.04 && pos_y < 0.16   
                Acalc = Aij_front_vert(A, j, i, miferro, 0, mibobina, Jz(pos_y));

            % Fronteira vertical ar-bobina
            elseif pos_x == 0.14 && pos_y > 0.04 && pos_y < 0.16 
                Acalc = Aij_front_vert(A, j, i, miar, 0, mibobina, Jz(pos_y));
            
            % Fronteira vertical ar-ferro esquerda
            elseif pos_x == 0.04
                Acalc = Aij_front_vert(A, j, i, miferro, 0, miar, 0);
            
            % Fronteira vertical ar-ferro direita
            elseif (pos_y < 0.04 || pos_y > 0.16) && pos_x == 0.05
                Acalc = Aij_front_vert(A, j, i, miar, 0, miferro, 0);
            
            % Fronteira horizontal bobina-ferro cima
            elseif pos_y == 0.16 && pos_x > 0.14 && pos_x < 0.16
                Acalc = Aij_front_hori(A, j, i, mibobina, Jz(pos_y), miferro, 0);
            
            % Fronteira horizontal bobina-ferro baixo
            elseif pos_y == 0.04 && pos_x > 0.14 && pos_x < 0.16
                Acalc = Aij_front_hori(A, j, i, miferro, 0, mibobina, Jz(pos_y));

            % Fronteira horizontal ar-ferro cima
            elseif pos_y == 0.16 && pos_x > 0.05 && pos_x < 0.14
                Acalc = Aij_front_hori(A, j, i, miar, 0, miferro, 0);
            
            % Fronteira horizontal ar-ferro baixo
            elseif pos_y == 0.04 && pos_x > 0.05 && pos_x < 0.14
                Acalc = Aij_front_hori(A, j, i, miferro, 0, miar, 0);
            
            % Quinas bobina-ferro-ar ???
            
            % Interior do domínio na bobina
            elseif pos_y > 0.04 && pos_y < 0.16 && ((pos_x > 0.14 && pos_x < 0.16) || pos_x > 0.20)      
                Acalc = Aij_interior(A, j, i, mibobina, Jz(pos_y));
                
            % Interior do domínio ferro ou ar                      
            else
                Acalc = Aij_interior(A, j, i, 0, 0);
                
            end
            
            % Calcula valor novo de A no ponto por Liebmann
            A(j, i) = lambda*Acalc + (1-lambda)*Aold(j, i);

            % Calcula erro no ponto e substitui o erro máximo se for maior
            erro = abs((A(j, i) - Aold(j, i))/A(j, i));
            if erro > erro_max
                erro_max = erro;
            end
                  
        end
    end
    
end

%[Y,X] = meshgrid(0:dx:0.21,0:dx:0.23);