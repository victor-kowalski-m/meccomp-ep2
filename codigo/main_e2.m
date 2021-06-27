clear
clc
close all

%% ESCOLHA DO ITEM
item = "e2"; % "ad", "e1" ou "e2"

%% Parâmetros iniciais
dx = 0.0025;
dt_min = dx^2;
[dx, dy, rows, cols, row_eq, col_eq, A, JZ, Fronteiras, vertical, ...
    horizontal, vazio_direita, mi0, MIx, MIy, Sigma, dt, tempos]...
     = gera_parametros_iniciais(item, dx, dt_min);

%% Funções auxiliares para cálculo de Aij
% equacoes = gera_equacoes(item, dx, dt);

equacoes = cell(0);
dx2 = dx^2;


% Em pontos da bobina %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aij em um ponto no interior da bobina
equacoes{end + 1} = @(A, j, i, k, t, mix, miy, Jz, sigma)...
    (dx.^(-2).*mix.^(-1).*miy.^(-1).*sigma.^(-1).*(dt.* ...
    miy.*(A((-1)+j,i,k)+(-2).*A(j,i,k)+A(1+j,i,k))+mix.*(dt.*(A(j,(-1) ...
    +i,k)+(-2).*A(j,i,k)+A(j,1+i,k))+dx2.*miy.*(sigma.*A(j,i,k)+dt.* ...
    Jz.*cos(60.*t)))));

% Aij na fronteira vertical da bobina
equacoes{end + 1} = @(A, j, i, k, t, mix1, miy1, Jz1, mix2, ...
    miy2, Jz2, sigma1, sigma2)...
    (dx.^(-2).*mix1.^(-1).*mix2.^(-1).*miy1.^(-1).* ...
    miy2.^(-1).*(sigma1+sigma2).^(-1).*(((-2).*dt.*mix2.*miy1.*miy2+ ...
    mix1.*((-2).*dt.*miy1.*miy2+mix2.*((-2).*dt.*miy2+miy1.*((-2).*dt+ ...
    dx2.*miy2.*(sigma1+sigma2))))).*A(j,i,k)+dt.*(2.*mix1.*mix2.* ...
    miy1.*A(j,1+i,k)+miy2.*(2.*mix1.*mix2.*A(j,(-1)+i,k)+miy1.*((mix1+ ...
    mix2).*A((-1)+j,i,k)+(mix1+mix2).*A(1+j,i,k)+dx2.*(Jz1+Jz2).* ...
    mix1.*mix2.*cos(60.*t))))));

% Aij na fronteira horizontal da bobina
equacoes{end + 1} = @(A, j, i, k, t, mix1, miy1, Jz1, mix2, ...
    miy2, Jz2, sigma1, sigma2) ...
    (dx.^(-2).*mix1.^(-1).*mix2.^(-1).*miy1.^(-1).* ...
    miy2.^(-1).*(sigma1+sigma2).^(-1).*(dt.*mix1.*mix2.*(miy1+miy2).* ...
    A(j,(-1)+i,k)+((-2).*dt.*mix2.*miy1.*miy2+mix1.*((-2).*dt.*miy1.* ...
    miy2+mix2.*((-2).*dt.*miy2+miy1.*((-2).*dt+dx2.*miy2.*(sigma1+ ...
    sigma2))))).*A(j,i,k)+dt.*(mix1.*mix2.*(miy1+miy2).*A(j,1+i,k)+ ...
    miy1.*miy2.*(2.*mix1.*A(1+j,i,k)+mix2.*(2.*A((-1)+j,i,k)+dx2.*( ...
    Jz1+Jz2).*mix1.*cos(60.*t))))));


% Aij um ponto no interior do domínio
equacoes{end + 1} = @(A, j, i, k, mix, miy, Jz)...
    ((1/2).*(mix+miy).^(-1).*(dx2.*Jz.*mix.*miy+mix.*...
    (A(j,i-1, k+1)+A(j,i+1, k+1))+miy.*(A(j-1,i, k+1)+A(j+1,i, k+1))));

% Aij na fronteira vertical entre dois meios
equacoes{end + 1} = @(A, j, i, k, mix1, miy1, Jz1, mix2, miy2, Jz2)...
    ((1/2).*(mix2.*miy1.*miy2+mix1.*(miy1.*miy2+mix2.*...
    (miy1+miy2))).^(-1).*(2.*mix1.*mix2.*miy1.*A(j,i+1, k+1)+miy2.*...
    (2.*mix1.*mix2.*A(j,i-1, k+1)+miy1.*(dx2.*(Jz1+Jz2).*...
    mix1.*mix2+(mix1+mix2).*A(j-1,i, k+1)+(mix1+mix2).*A(j+1,i, k+1)))));

% Aij na fronteira horizontal entre dois meios
equacoes{end + 1} = @(A, j, i, k, mix1, miy1, Jz1, mix2, miy2, Jz2) ...
    ((1/2).*(mix2.*miy1.*miy2+mix1.*(miy1.*miy2+mix2.*...
    (miy1+miy2))).^(-1).*(mix1.*mix2.*(miy1+miy2).*A(j,i-1, k+1)+...
    mix1.*mix2.*(miy1+miy2).*A(j,i+1, k+1)+miy1.*miy2.*...
    (mix2.*(dx2.*(Jz1+Jz2).*mix1+2.*A(j-1,i, k+1))+2.*mix1.*A(j+1,i, k+1))));

%% Aplicação do MDF
% [A, iters] = ...
%     aplica_MDF(item, A, rows, cols, row_eq, col_eq, Fronteiras, ...
%     vertical, horizontal, vazio_direita, equacoes, MIx, MIy, JZ, Sigma, dt, tempos);

% Equações na bobina
Aij_int_bobina = equacoes{1};
Aij_vert_bobina = equacoes{2};
Aij_hori_bobina = equacoes{3};

% Equações fora da bobina
Aij_int_fora = equacoes{4};
Aij_vert_fora = equacoes{5};
Aij_hori_fora = equacoes{6};

instantes = tempos(end)-1;
Tempo = 0:dt:instantes*dt; % Vetor de tempos
A = repmat(A,[1 1 instantes+1]); % Matriz inicial de valores de A no tempo

rows_bobina = row_eq(4)+1:row_eq(16)-1;
cols_bobina = [col_eq(14)+1:col_eq(16)-1 col_eq(20)+1:col_eq(22)-1];

iters = 0; % contador de iteracoes

for k = 1:instantes
            
    iters = iters + 1;
    A(:, :, k+1) = A(:, :, k);
    Ak = A(:, :, k+1);
    asda = 0;

    % Itera bobina por linha e coluna
    for j=rows_bobina 
        for i=cols_bobina

            asda = 0;

            % Se está embaixo ou em cima da bobina externa, pula linha
            if Fronteiras(j, i) == vazio_direita
                break;
            end

            % Escolhe qual equação de Aij utilizar
            if Fronteiras(j, i) == vertical % se é fronteira vertical
                A(j, i, k+1) = Aij_vert_bobina( ...
                    A, j, i, k, Tempo(k), ...
                    MIx(j, i-1), MIy(j, i-1), JZ(j,i-1), ...
                    MIx(j, i+1), MIy(j, i+1), JZ(j,i+1), ...
                    Sigma(j, i-1), Sigma(j, i+1));

            elseif Fronteiras(j, i) == horizontal % se é fronteira horizontal
                A(j, i, k+1) = Aij_hori_bobina( ...
                    A, j, i, k, Tempo(k), ...
                    MIx(j-1, i), MIy(j-1, i), JZ(j-1,i), ...
                    MIx(j+1, i), MIy(j+1, i), JZ(j+1,i), ...
                    Sigma(j-1, i), Sigma(j+1, i));

            else % se estiver em um ponto interior do domínio
                A(j, i, k+1) = Aij_int_bobina( ...
                    A, j, i, k, Tempo(k), ...
                    MIx(j, i), MIy(j, i), JZ(j, i), Sigma(j, i));
            end                                    

        end
    end
    
    
    Aold = A; % salva estado da matriz A antes de fazer o novo cálculo
    lambda = 1.75;
    for j=2:rows-1
        for i=2:cols-1

            % Se está embaixo ou em cima da bobina externa, pula linha
            if (Fronteiras(j, i) == vazio_direita)
                break;
            end
            
            if ((j > row_eq(4) && j < row_eq(16)) && ...
                ((i > col_eq(14) && i < col_eq(16)) || ...
                (i > col_eq(20) && i < col_eq(22))))
                continue;
            end
                        
            % Escolhe qual equação de Aij utilizar
            if Fronteiras(j, i) == vertical % se é fronteira vertical
                Acalc = Aij_vert_fora(A, j, i, k, ...
                MIx(j,i-1), MIy(j,i-1), JZ(j,i-1), ...
                MIx(j,i+1), MIy(j,i+1), JZ(j,i+1));
            elseif Fronteiras(j, i) == horizontal % se é fronteira horizontal
                Acalc = Aij_hori_fora(A, j, i, k, ...
                MIx(j-1,i), MIy(j-1,i), JZ(j-1,i), ...
                MIx(j+1,i), MIy(j+1,i), JZ(j+1,i));
            else % se estiver em um ponto interior do domínio
                Acalc = Aij_int_fora(A, j, i, k, ...
                MIx(j, i), MIy(j, i), JZ(j, i));
            end                                   

            % Calcula valor novo de A no ponto por Liebmann
            A(j, i, k+1) = lambda*Acalc + (1-lambda)*Aold(j, i, k);
            
            if i == col_eq(18)
                asda = A(j, i, k+1);
                xasd=0;
            end

        end
    end
end

%% Cálculo do vetor B e H
[Bx, By, Hx, Hy, Fela_x, Fela_y]...
    = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, cols, ...
    col_eq, Fronteiras, vazio_direita, tempos);

%% Plots
plota(A, dx, dy, Bx, By, Hx, Hy, tempos);

