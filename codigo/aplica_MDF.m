function [A, iter, Bx, By, Hx, Hy, Fela_x, Fela_y] = ...
    aplica_MDF(item, A, rows, cols, row_eq, col_eq, Fronteiras, ...
    vertical, horizontal, vazio_direita, equacoes, MIx, MIy, JZ, Sigma,...
    dt, tempos, mi0, dx, dy)
    
    % Se é o item a-d ou e1, resolve por Liebmann
    if any(item == ["ad" "e1"])
        
        % Nomeia as equações recebidas como argumento
        Aij_interior = equacoes{1};
        Aij_vert = equacoes{2};
        Aij_hori = equacoes{3};
        
        lambda = 1.75; % lambda de sobrerrelaxação
        tolerancia = 0.0001; % criterio de parada do erro
        limite_iters = 100000; % numero maximo de iterações
        iter = 0; % contador de iteracoes
        erro_max = inf; % erro inicial para entrar no loop

        % Calcula A por liebmann até alcançar a tolerância desejada ou 
        % bater o limite de iterações
        while erro_max > tolerancia && iter < limite_iters

            iter = iter + 1; % incrementa contador
            erro_max = 0; % zera erro_maximo
            Aold = A; % salva estado da matriz A antes de fazer o novo cálculo
            
            % Itera matriz por linha e coluna
            for j=2:rows-1
                for i=2:cols-1

                    % Se está embaixo ou em cima da bobina externa, pula linha
                    if Fronteiras(j, i) == vazio_direita
                        break;
                    end

                    % Escolhe qual equação utilizar para calcular Aij
                    if Fronteiras(j, i) == vertical % se é fronteira vertical
                        Acalc = Aij_vert(A, j, i, ...
                        MIx(j,i-1), MIy(j,i-1), JZ(j,i-1), ...
                        MIx(j,i+1), MIy(j,i+1), JZ(j,i+1));
                    elseif Fronteiras(j, i) == horizontal % se é fronteira horizontal
                        Acalc = Aij_hori(A, j, i, ...
                        MIx(j-1,i), MIy(j-1,i), JZ(j-1,i), ...
                        MIx(j+1,i), MIy(j+1,i), JZ(j+1,i));
                    else % se estiver em um ponto interior do domínio
                        Acalc = Aij_interior(A, j, i, ...
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
        
        % Calcula grandezas B, H e Fela com a matriz A encontrada
        [Bx, By, Hx, Hy, Fela_x, Fela_y]...
            = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, cols, ...
            col_eq, Fronteiras, vazio_direita, 1, 1);
        
    
    % Se é o item e2, resolve usando método explícito na bobina 
    elseif item == "e2"
        
        % Nomeia as equações recebidas como argumento
        % Equações na bobina
        Aij_int_bobina = equacoes{1};
        Aij_vert_bobina = equacoes{2};
        Aij_hori_bobina = equacoes{3};
        % Equações fora da bobina
        Aij_int_fora = equacoes{4};
        Aij_vert_fora = equacoes{5};
        Aij_hori_fora = equacoes{6};
                        
        % Pontos dentro da bobina
        rows_bobina = row_eq(4):row_eq(16);
        cols_bobina = [col_eq(14):col_eq(16) col_eq(20):col_eq(22)-1];
                
        Tempo = 0:dt:tempos(end)*dt; % Vetor de tempos
        iter = 0; % contador de iteracoes
        
        Atempos = repmat(A, [1 1 length(tempos)]); % Matriz que guardará A nos tempos desejados
        A = repmat(A, [1 1 2]); % Guarda o estado atual e anterior da matriz A
        
        % Matrizes que guardarão B e H nos tempos desejados
        [Bx, By, Hx, Hy] = deal(zeros(rows, cols, length(tempos)));
        
        % Vetores que guardarão Fela para todos os tempos
        [Fela_x, Fela_y] = deal(zeros(1, tempos(end)));
        
        k = 1; % indice temporal de A, k antigo e k+1 atual 
        lambda = 1.75; % lambda de sobrerrelaxação
        
        % Itera para cada tempo até o último tempo desejado
        for t=Tempo
            
            A(: , :, k) = A(:, :, k+1); % Matriz atual da iteração anterior passa a ser a antiga
            iter = iter + 1; % incrementa contador de iterações
            
            % Itera bobina por linha e coluna
            for j=rows_bobina 
                for i=cols_bobina
                    
                    % Se está embaixo ou em cima da bobina externa, pula linha
                    if Fronteiras(j, i) == vazio_direita
                        break;
                    end
                    
                    % Se é um dos cantos da bobina, pula
                    if any([row_eq(4) row_eq(16)] == j) && ...
                            any([col_eq(14) col_eq(16) col_eq(20)] == i)
                        continue
                    end
                    
                    % Escolhe qual equação de Aij (na bobina) utilizar
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
                    else % se é um ponto interior
                        A(j, i, k+1) = Aij_int_bobina( ...
                            A, j, i, k, t, ...
                            MIx(j, i), MIy(j, i), JZ(j, i), Sigma(j, i));                                   
                    end
                        
                end
            end
                      
            % Itera pontos fora da bobina
            for j=2:rows-1
                for i=2:col_eq(20)-1

                    % Se for um ponto da bobina pula (exceto cantos)
                    if i >= col_eq(14) && i <= col_eq(16) && ...
                        j >= row_eq(4) && j <= row_eq(16) && ...
                        ~(any([row_eq(4) row_eq(16)] == j) && ...
                            any([col_eq(14) col_eq(16) col_eq(20)] == i))
                        continue
                    end
                                       
                    % Escolhe qual equação de Aij utilizar
                    if Fronteiras(j, i) == vertical % se é fronteira vertical
                        Acalc = Aij_vert_fora(A, j, i, k, t,...
                        MIx(j,i-1), MIy(j,i-1), JZ(j,i-1), ...
                        MIx(j,i+1), MIy(j,i+1), JZ(j,i+1));
                    elseif Fronteiras(j, i) == horizontal % se é fronteira horizontal
                        Acalc = Aij_hori_fora(A, j, i, k, t,...
                        MIx(j-1,i), MIy(j-1,i), JZ(j-1,i), ...
                        MIx(j+1,i), MIy(j+1,i), JZ(j+1,i));
                    else % se estiver em um ponto interior do domínio
                        Acalc = Aij_int_fora(A, j, i, k, t,...
                        MIx(j, i), MIy(j, i), JZ(j, i));
                    end                                   

                    % Calcula valor novo de A no ponto com sobrerrelaxação
                    A(j, i, k+1) = lambda*Acalc + (1-lambda)*A(j, i, k);

                end
            end
           
            % Calcula B, H e Fela para o instante de tempo
            [Bx_t, By_t, Hx_t, Hy_t, Fela_x_t, Fela_y_t]...
                = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, ...
                cols, col_eq, Fronteiras, vazio_direita, k+1, 0);
            
            % Adiciona Fela na array
            Fela_x(iter) =  Fela_x_t;
            Fela_y(iter) =  Fela_y_t;
            
            % Se for um dos tempos desejados, salva A, B e H
            for idx_t=1:length(tempos)
                if iter == tempos(idx_t)
                    Atempos(:, :, idx_t) = A(:, :, k+1);
                    Bx(:, :, idx_t)= Bx_t;
                    By(:, :, idx_t)= By_t;
                    Hx(:, :, idx_t)= Hx_t;
                    Hy(:, :, idx_t)= Hy_t;
                end
            end
            
        end
        
        A = Atempos; % renomeia Atempos para A
        
    end  
end

