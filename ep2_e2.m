clear
clc
close all

%% Parâmetros iniciais para o item e)
[dx, dy] = deal(0.0025); % passos em metros

rows = 0.20/dy + 1; % numero de linhas da matriz
cols = 0.22/dx + 1; % numero de colunas da matriz
A0 = zeros(rows, cols); % matriz inicial com zeros

% Funções para calcular os pontos da malha equivalentes a x e y no desenho
col_eq = @(x) x/22*(cols-1)+1; 
row_eq = @(y) y/20*(rows-1)+1;

% Pontos de início de fim da bobina externa
ini_bobina = [row_eq(4) col_eq(20)]; % linha e coluna iniciais da bobina
fim_bobina = [row_eq(16) col_eq(22)]; % linha e coluna finais da bobina

% Define matriz de permeabilidade magnética no domínio
mi0 = 4*pi*10^-7; % permeabilidade magnética do vacuo, do ar e da bobina
miferrox = 1200*mi0; % permeabilidade magnética do ferro em x
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
    for i = col_eq(20)+1:col_eq(22)
        JZ(j, i) = Jz((j-1)*dy);
    end
    for i = col_eq(14):col_eq(16)-1 
        JZ(j, i) = -Jz((j-1)*dy);  
    end
end

%% Funções auxiliares para cálculo de Aij

dx2 = dx^2;
sigma = 4*10^6;
dt = 0.00001; % Passo no tempo
lenT = 100; % Número de passos
T = 0:dt:(lenT-1)*dt; % Vetor de tempos
A = repmat(A0,[1 1 lenT]); % Matriz inicial de valores de A no tempo

% Função de cálculo de Aij com mi diferente em x e y em um ponto no interior do domínio
Aijkmais1_interior = @(A, j, i, k, t, mix, miy, Jz)...
    (dx.^(-2).*mix.^(-1).*miy.^(-1).*sigma.^(-1).*(dt.*miy.* ...
    (A(k,j-1,i)+(-2).*A(k,j,i)+A(k,j+1,i))+mix.* ...
    (dt.*(A(k,j,i-1)+(-2).*A(k,j,i)+A(k,j,i+1))+dx2.*miy.* ...
    (sigma.*A(k,j,i)+dt.*Jz.*cos(60.*t)))));

% Função de cálculo de Aij com mi diferente em x e y na fronteira vertical entre dois meios
Aijkmais1_front_vert = @(A, j, i, k, t, mix1, miy1, Jz1, mix2, miy2, Jz2)...
    ((1/2).*dx.^(-2).*mix1.^(-1).*mix2.^(-1).*miy1.^(-1).* ...
    miy2.^(-1).*sigma.^(-1).*((-2).*(dt.*mix2.*miy1.*miy2+mix1.* ...
    (dt.*miy1.*miy2+mix2.*(dt.*miy2+miy1.*(dt+(-1).*dx2.*miy2.*sigma)))).* ...
    A(k,j,i)+dt.*(2.*mix1.*mix2.*miy1.*A(k,j,i+1)+miy2.* ...
    (2.*mix1.*mix2.*A(k,j,i-1)+miy1.*((mix1+mix2).* ...
    A(k,j-1,i)+(mix1+mix2).*A(k,j+1,i)+dx2.* ...
    (Jz1+Jz2).*mix1.*mix2.*cos(60.*t))))));

% Função de cálculo de Aij com mi diferente em x e y na fronteira horizontal entre dois meios
Aijkmais1_front_hori = @(A, j, i, k, t, mix1, miy1, Jz1, mix2, miy2, Jz2) ...
    ((1/2).*dx.^(-2).*mix1.^(-1).*mix2.^(-1).*miy1.^(-1).* ...
    miy2.^(-1).*sigma.^(-1).*(dt.*mix1.*mix2.*(miy1+miy2).* ...
    A(k,j,i-1)+(-2).*(dt.*mix2.*miy1.*miy2+mix1.* ...
    (dt.*miy1.*miy2+mix2.*(dt.*miy2+miy1.* ...
    (dt+(-1).*dx2.*miy2.*sigma)))).*A(k,j,i)+dt.* ...
    (mix1.*mix2.*(miy1+miy2).*A(k,j,i+1)+miy1.*miy2.* ...
    (2.*mix1.*A(k,j+1,i)+mix2.*(2.*A(k,j-1,i)+dx2.* ...
    (Jz1+Jz2).*mix1.*cos(60.*t))))));
 
%% Aplicação do MDF

for k = 1:lenT - 1
        
    % Itera matriz por linha e coluna
    for j=2:rows-1
        for i=2:cols-1
                       
            % Se está embaixo ou em cima da bobina externa, pula
            if i >= ini_bobina(2) && (j <= ini_bobina(1) || j >= fim_bobina(1))
                continue;
            end
            
            % Escolhe qual equação de Aij utilizar
            if MIx(j, i-1) ~= MIx(j, i+1) % se é fronteira vertical
                A(k+1, j, i) = Aijkmais1_front_vert( ...
                    A, j, i, k, T(k), ...
                    MIx(j, i-1), MIy(j, i-1), JZ(j,i-1), ...
                    MIx(j, i+1), MIy(j, i+1), JZ(j,i+1));
                
            elseif MIx(j-1,i) ~= MIx(j+1,i) % se é fronteira horizontal
                A(k+1, j, i) = Aijkmais1_front_hori( ...
                    A, j, i, k, T(k), ...
                    MIx(j, i-1), MIy(j, i-1), JZ(j,i-1), ...
                    MIx(j, i+1), MIy(j, i+1), JZ(j,i+1));
                
            else % se estiver em um ponto interior do domínio
                A(k+1, j, i) = Aijkmais1_interior( ...
                    A, j, i, k, T(k), ...
                    MIx(j, i), MIy(j, i), JZ(j, i));
                    
            end                                    
                              
        end
    end
    
end

% %% Cálculo do vetor B e H
% 
% % Inicialização das matrizes de densidade superficial de fluxo magnético e
% % intensidade de campo magnético
% [Bx, By, Hx, Hy] = deal(zeros(rows, cols));
% 
% for j = 1:rows
%    for i = 1:cols
%        
%        pos_x = (i-1)*dx;
%        pos_y = (j-1)*dy;
%        
%        % Se for canto ele pula
%        if ((pos_y == 0 && pos_x == 0) || (pos_y == 0.20 && pos_x == 0) || ...
%                (pos_y == 0 && pos_x == 0.22) || (pos_y == 0.20 && pos_x == 0.22))
%            continue;
%        end
%        
%        % Borda inferior
%        if (pos_y == 0 && pos_x > 0 && pos_x < 20) ||...
%                (pos_y == 0.04 && pos_x > 0.20 && pos_x < 0.22)
%            Bx(j,i) = (-A(j+2,i) + 4*A(j+1,i) - 3*A(j,i))/(2*dy);
%            
%        % Borda superior
%        elseif (pos_y == 0.20 && pos_x > 0 && pos_x < 20) ||...
%                (pos_y == 0.16 && pos_x > 0.20 && pos_x < 0.22)
%            Bx(j,i) = (3*A(j,i) - 4*A(j-1,i) + A(j-2,i))/(2*dy);
%            
%        % Borda esquerda
%        elseif pos_x == 0   
%            By(j,i) = -(-A(j,i+2) + 4*A(j,i+1) - 3*A(j,i))/(2*dx);
%            
%        % Borda direita
%        elseif pos_x == 0.22
%            By(j,i) = -(3*A(j,i) - 4*A(j,i-1) + A(j,i-2))/(2*dx);
%            
%        % Parte interna
%        else
%            Bx(j,i) = (A(j+1,i) - A(j-1,i))/(2*dy);
%            By(j,i) = -(A(j,i+1) - A(j,i-1))/(2*dy);
%        end
%        
%        Hx(j,i) = Bx(j,i)/MIx(j,i);
%        Hy(j,i) = By(j,i)/MIy(j,i);
%        
%    end
% end

%% Plots

% Plot de Az
[X,Y] = meshgrid(0:dx:0.22,0:dy:0.20);
surf(X,Y,A, 'LineStyle', ':')
colorbar

% % Plot de B
% figB = figure("Name", "Vetor de densidade superficial de fluxo magnético");
% fluxo = quiver(X,Y,Bx,By);
% axis equal
% xlabel("x(m)")
% ylabel("y(m)")
% grid()
% 
% % Plot de H
% figH = figure("Name", "Vetor de intensidade de campo magnético");
% campo = quiver(X,Y,Hx,Hy);
% axis equal
% xlabel("x(m)")
% ylabel("y(m)")
% grid()
% 
% %% Força
% 
% Fela_x = 1/(2*mi0)*(...
%     trapz(0:dx:0.2, Bx(:, col_eq(4)).^2 - By(:, col_eq(4)).^2) -...
%     trapz(0:dx:0.04, 2*Bx(row_eq(20), col_eq(0):col_eq(4)) ...
%     .*By(row_eq(20), col_eq(0):col_eq(4))) +...
%     trapz(0:dx:0.2, Bx(:, col_eq(0)).^2 - By(:, col_eq(0)).^2) -...
%     trapz(0:dx:0.04, 2*Bx(row_eq(0), col_eq(0):col_eq(4)) ...
%     .*By(row_eq(0), col_eq(0):col_eq(4))) );
% 
% Fela_y = 1/(2*mi0)*(...
%     trapz(0:dx:0.2, 2*Bx(:, col_eq(4)).*By(:, col_eq(4))) -...
%     trapz(0:dx:0.04, By(row_eq(20), col_eq(0):col_eq(4)).^2 ...
%     - Bx(row_eq(20), col_eq(0):col_eq(4)).^2) +...
%     trapz(0:dx:0.2, 2*Bx(:, col_eq(0)).*By(:, col_eq(0))) -...
%     trapz(0:dx:0.04, By(row_eq(0), col_eq(0):col_eq(4)).^2 ...
%     - Bx(row_eq(0), col_eq(0):col_eq(4)).^2));






