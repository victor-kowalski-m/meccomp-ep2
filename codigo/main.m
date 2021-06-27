clear
clc
close all

%% ESCOLHA DO ITEM
item = "e2"; % "ad", "e1" ou "e2"

%% Parâmetros iniciais
dx = 0.001;
dt_min = dx^2;
[dx, dy, rows, cols, row_eq, col_eq, A, JZ, Fronteiras, vertical, ...
    horizontal, vazio_direita, mi0, MIx, MIy, Sigma, dt, tempos]...
     = gera_parametros_iniciais(item, dx, dt_min);

%% Funções auxiliares para cálculo de Aij
equacoes = gera_equacoes(item, dx, dt);

%% Aplicação do MDF
[A, iters] = ...
    aplica_MDF(item, A, rows, cols, row_eq, col_eq, Fronteiras, ...
    vertical, horizontal, vazio_direita, equacoes, MIx, MIy, JZ, Sigma, dt, tempos);

%% Cálculo do vetor B e H
[Bx, By, Hx, Hy, Fela_x, Fela_y]...
    = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, cols, ...
    col_eq, Fronteiras, vazio_direita, tempos);

%% Plots
plota(A, dx, dy, Bx, By, Hx, Hy, tempos);

