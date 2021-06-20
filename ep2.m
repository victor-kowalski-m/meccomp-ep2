clear
clc
close all

%% Parâmetros iniciais
[dx, dy] = deal(0.001); % passos em metros

rows = 0.20/dy + 1; % numero de linhas da matriz
cols = 0.22/dx + 1; % numero de colunas da matriz
A = zeros(rows, cols); % matriz inicial com zeros

% Funções para calcular os pontos da malha equivalentes a x e y no desenho
col_eq = @(x) x/22*(cols-1)+1; 
row_eq = @(y) y/20*(rows-1)+1;

% Pontos de início de fim da bobina externa
ini_bobina = [row_eq(4) col_eq(20)]; % linha e coluna iniciais da bobina
fim_bobina = [row_eq(16) col_eq(22)]; % linha e coluna finais da bobina

% Define matriz de permeabilidade magnética no domínio
mi0 = 4*pi*10^-7; % permeabilidade magnética do vacuo, do ar e da bobina
miferro = 2500*mi0; % permeabilidade magnética do ferro
MI = ones(rows, cols)*mi0;
MI(:, col_eq(0):col_eq(4)) = miferro;
MI(:, col_eq(16):col_eq(20)) = miferro;
MI(row_eq(0):row_eq(4), col_eq(5):col_eq(20)) = miferro;
MI(row_eq(16):row_eq(20), col_eq(5):col_eq(20)) = miferro;
MI(row_eq(4), col_eq(5)) = mi0;
MI(row_eq(16), col_eq(5)) = mi0;

% Define matriz de e densidade superficial de corrente elétrica no domínio
JZ = zeros(rows, cols);
Jz = @(y) (2*10^6*cos(pi*y/(12*10^-2)) + 8*10^5);
for j=row_eq(4)+1:row_eq(16)-1 
    for i=cat(2, col_eq(14):col_eq(16)-1, col_eq(20)+1:col_eq(22))    
        JZ(j, i) = Jz((j-1)*dy);
    end
end

%% Funções auxiliares para cálculo de Aij

% Função de cálculo de Aij em um ponto no interior do domínio
Aij_interior = @(A, j, i, mi, Jz)...
    (A(j, i+1) + A(j, i-1) + A(j+1, i) + A(j - 1, i))/4 + mi*Jz;

% Função de cálculo de Aij na fronteira vertical entre dois meios
Aij_front_vert = @(A, j, i, mi1, Jz1, mi2, Jz2)...
    (2*(mi2*A(j,i-1) + mi1*A(j,i+1)) + (mi1 + mi2)*(A(j-1,i) + A(j+1,i))...
    + (mi1*mi2*dx^2)*(Jz1 + Jz2))/(4*(mi1 + mi2)); 

% Função de cálculo de Aij na fronteira horizontal entre dois meios
Aij_front_hori = @(A, j, i, mi1, Jz1, mi2, Jz2) ...
    (2*(mi2*A(j-1,i) + mi1*A(j+1,i)) + (mi1 + mi2)*(A(j,i-1) + A(j,i+1))...
    + (mi1*mi2*dx^2)*(Jz1 + Jz2))/(4*(mi1 + mi2)); 
 
%% Aplicação do MDF

lambda = 1.75; % lambda de sobrerrelaxação
tolerancia = 0.0001; % criterio de parada do erro
limite_iters = 100000; % numero maximo de iterações
iters = 0; % contador de iteracoes
erro_max = inf; % erro inicial para entrar no loop

% Calcula A por liebmann até alcançar a tolerância desejada ou bater o
% limite de iterações
while erro_max > tolerancia && iters < limite_iters
    
    iters = iters + 1; % incrementa contador
    erro_max = 0; % zera erro_maximo
    Aold = A; % salva estado da matriz A antes de fazer o novo cálculo
    
    % Itera matriz por linha e coluna
    for j=2:rows-1
        for i=2:cols-1
                       
            % Se está embaixo ou em cima da bobina externa, pula
            if i >= ini_bobina(2) && (j <= ini_bobina(1) || j >= fim_bobina(1))
                continue;
            end
            
            % Escolhe qual equação de Aij utilizar
            if MI(j, i-1) ~= MI(j, i+1) % se é fronteira vertical
                Acalc = Aij_front_vert(A, j, i, ...
                    MI(j,i-1), JZ(j,i-1), MI(j,i+1), JZ(j,i+1));
            elseif MI(j-1,i) ~= MI(j+1,i) % se é fronteira horizontal
                Acalc = Aij_front_hori(A, j, i, ...
                    MI(j-1,i), JZ(j-1,i), MI(j+1,i), JZ(j+1,i));
            else % se estiver em um ponto interior do domínio
                Acalc = Aij_interior(A, j, i, MI(j, i), JZ(j, i));
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

%% Plots

[X,Y] = meshgrid(0:dx:0.22,0:dy:0.20);
surf(X,Y,A)
colorbar