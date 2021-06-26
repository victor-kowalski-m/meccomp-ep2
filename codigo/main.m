clear
clc
close all

%% ESCOLHA DO ITEM
item = "e1"; % "ad", "e1" ou "e2"
lenT = 0; % Número de passos no tempo para item e2

%% Parâmetros iniciais
[dx, dy, rows, cols, col_eq, row_eq, A, JZ, Fronteiras, ...
     vertical, horizontal, vazio_direita, mi0, MIx, MIy, Sigma]...
     = gera_parametros_iniciais(item, 0.0025);

%% Funções auxiliares para cálculo de Aij
[Aij_interior, Aij_front_vert, Aij_front_hori] = ...
    gera_equacoes(item, dx);

%% Aplicação do MDF
[A, iters, erro_max] = ...
    aplica_MDF(item, A, rows, cols, Fronteiras, vertical, horizontal, ...
    vazio_direita, Aij_interior, Aij_front_vert, Aij_front_hori, ...
    MIx, MIy, JZ, lenT);

%% Cálculo do vetor B e H
[Bx, By, Hx, Hy, Fela_x, Fela_y]...
    = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, cols, ...
    col_eq, Fronteiras, vazio_direita);

%% Plots
[figA, figB, figH] = plota(A, dx, dy, Bx, By, Hx, Hy);

