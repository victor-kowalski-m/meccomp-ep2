clear
clc
close all

%% Parâmetros iniciais para o item e)
[dx, dy] = deal(0.001); % passos em metros

rows = 0.20/dy + 1; % numero de linhas da matriz
cols = 0.22/dx + 1; % numero de colunas da matriz
A = zeros(rows, cols); % matriz inicial com zeros

% Funções para calcular os pontos da malha equivalentes a x e y no desenho
col_eq = @(x) x/22*(cols-1)+1; 
row_eq = @(y) y/20*(rows-1)+1;

% Define matriz de permeabilidade magnética no domínio
mi0 = 4*pi*10^-7; % permeabilidade magnética do vacuo, do ar e da bobina
miferrox = 1200*mi0; % 1200 permeabilidade magnética do ferro em x
miferroy = 2500*mi0; % permeabilidade magnética do ferro em y
[MIx,MIy] = deal(ones(rows, cols)*mi0);
MIx(:, col_eq(0):col_eq(4)) = miferrox;
MIx(:, col_eq(16):col_eq(20)) = miferrox;
MIx(row_eq(0):row_eq(4), col_eq(5):col_eq(20)) = miferrox;
MIx(row_eq(16):row_eq(20), col_eq(5):col_eq(20)) = miferrox;
MIy(:, col_eq(0):col_eq(4)) = miferroy;
MIy(:, col_eq(16):col_eq(20)) = miferroy;
MIy(row_eq(0):row_eq(4), col_eq(5):col_eq(20)) = miferroy;
MIy(row_eq(16):row_eq(20), col_eq(5):col_eq(20)) = miferroy;

% Define matriz de e densidade superficial de corrente elétrica no domínio
JZ = zeros(rows, cols);
Jz = @(y) (2*10^6*cos(pi*(y-0.1)/(12*10^-2)) + 8*10^5);
for j=row_eq(4)+1:row_eq(16)-1 
    for i = col_eq(20)+1:col_eq(22)-1
        JZ(j, i) = Jz((j-1)*dy);
    end
    for i = col_eq(14)+1:col_eq(16)-1 
        JZ(j, i) = -Jz((j-1)*dy);  
    end
end

% Define fronteiras do domínio: 0 interior, 1 vertical, 2 horizontal
vertical = 1;
horizontal = 2;
vazio_direita = 3;
Fronteiras = zeros(rows, cols);
Fronteiras(2:end-1, col_eq(4)) = vertical;
Fronteiras([2:row_eq(4) row_eq(16):end-1], col_eq(5)) = vertical;
Fronteiras(row_eq(4):row_eq(16), [col_eq(16) col_eq(20)]) = vertical;
Fronteiras([row_eq(4) row_eq(16)], col_eq(5):col_eq(16)) = horizontal;
Fronteiras([2:row_eq(4) row_eq(16):end-1], col_eq(20):col_eq(22)-1) = vazio_direita;

%% Funções auxiliares para cálculo de Aij

dx2 = dx^2;

% Função de cálculo de Aij com mi diferente em x e y em um ponto no interior do domínio
Aij_2mis_interior = @(A, j, i, mix, miy, Jz)...
    ((1/2).*(mix+miy).^(-1).*(dx2.*Jz.*mix.*miy+mix.*...
    (A(j,i-1)+A(j,i+1))+miy.*(A(j-1,i)+A(j+1,i))));

% Função de cálculo de Aij com mi diferente em x e y na fronteira vertical entre dois meios
Aij_2mis_front_vert = @(A, j, i, mix1, miy1, Jz1, mix2, miy2, Jz2)...
    ((1/2).*(mix2.*miy1.*miy2+mix1.*(miy1.*miy2+mix2.*...
    (miy1+miy2))).^(-1).*(2.*mix1.*mix2.*miy1.*A(j,i+1)+miy2.*...
    (2.*mix1.*mix2.*A(j,i-1)+miy1.*(dx2.*(Jz1+Jz2).*...
    mix1.*mix2+(mix1+mix2).*A(j-1,i)+(mix1+mix2).*A(j+1,i)))));

% Função de cálculo de Aij com mi diferente em x e y na fronteira horizontal entre dois meios
Aij_2mis_front_hori = @(A, j, i, mix1, miy1, Jz1, mix2, miy2, Jz2) ...
    ((1/2).*(mix2.*miy1.*miy2+mix1.*(miy1.*miy2+mix2.*...
    (miy1+miy2))).^(-1).*(mix1.*mix2.*(miy1+miy2).*A(j,i-1)+...
    mix1.*mix2.*(miy1+miy2).*A(j,i+1)+miy1.*miy2.*...
    (mix2.*(dx2.*(Jz1+Jz2).*mix1+2.*A(j-1,i))+2.*mix1.*A(j+1,i))));
 
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
                       
            % Se está embaixo ou em cima da bobina externa, pula linha
            if Fronteiras(j, i) == vazio_direita
                break;
            end
            
            % Escolhe qual equação de Aij utilizar
            if Fronteiras(j, i) == vertical % se é fronteira vertical
                Acalc = Aij_2mis_front_vert(A, j, i, ...
                    MIx(j,i-1), MIy(j,i-1), JZ(j,i-1), ...
                    MIx(j,i+1), MIy(j,i+1), JZ(j,i+1));
            elseif Fronteiras(j, i) == horizontal % se é fronteira horizontal
                Acalc = Aij_2mis_front_hori(A, j, i, ...
                    MIx(j-1,i), MIy(j-1,i), JZ(j-1,i), ...
                    MIx(j+1,i), MIy(j+1,i), JZ(j+1,i));
            else % se estiver em um ponto interior do domínio
                Acalc = Aij_2mis_interior(A, j, i, ...
                    MIx(j, i), MIy(j, i), JZ(j, i));
            end                                    
            
            % Calcula valor novo de A no ponto por Liebmann
            A(j, i) = lambda*Acalc + (1-lambda)*Aold(j, i);

            % Calcula erro no ponto e substitui o erro máximo se for maior
            if A(j, i) ~= 0  
                erro = abs((A(j, i) - Aold(j, i))/A(j, i));
                    if erro > erro_max
                        erro_max = erro;
                    end
            end
                  
        end
    end
    
end

%% Cálculo do vetor B e H

% Inicialização das matrizes de densidade superficial de fluxo magnético e
% intensidade de campo magnético
[Bx, By, Hx, Hy] = deal(zeros(rows, cols));

for j = 2:rows-1
   for i = 2:cols-1
       
        % Se está embaixo ou em cima da bobina externa, pula linha
        if Fronteiras(j, i) == vazio_direita
            break;
        end
           
       Bx(j,i) = (A(j+1,i) - A(j-1,i))/(2*dy);
       By(j,i) = -(A(j,i+1) - A(j,i-1))/(2*dx);
       Hx(j,i) = Bx(j,i)/MIx(j,i);
       Hy(j,i) = By(j,i)/MIy(j,i);
       
   end
end

%% Plots

% Plot de Az
[X,Y] = meshgrid(0:dx:0.22,0:dy:0.20);
surf(X,Y,A, 'LineStyle', ':')
colorbar

% Plot de B
figB = figure("Name", "Vetor de densidade superficial de fluxo magnético");
fluxo = quiver(X,Y,Bx,By);
axis equal
xlabel("x(m)")
ylabel("y(m)")
grid()

% Plot de H
figH = figure("Name", "Vetor de intensidade de campo magnético");
campo = quiver(X,Y,Hx,Hy);
axis equal
xlabel("x(m)")
ylabel("y(m)")
grid()

%% Força

Fela_x = 1/(2*mi0)*(...
    trapz(-0.1:dx:0.1, Bx(:, col_eq(4)).^2 - By(:, col_eq(4)).^2) );

Fela_y = 1/(2*mi0)*(...
    trapz(-0.1:dx:0.1, 2*Bx(:, col_eq(4)).*By(:, col_eq(4))) );


fprintf('A força Fx = %8.4f \n',Fela_x);
fprintf('A força Fy = %8.4f \n',Fela_y);



