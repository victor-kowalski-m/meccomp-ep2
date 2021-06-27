function [A, iter] = ...
    aplica_MDF(item, A, rows, cols, row_eq, col_eq, Fronteiras, ...
    vertical, horizontal, vazio_direita, equacoes, MIx, MIy, JZ, Sigma, dt, tempos)
    
    if any(item == ["ad" "e1"])
        
        Aij_interior = equacoes{1};
        Aij_vert = equacoes{2};
        Aij_hori = equacoes{3};
        
        lambda = 1.75; % lambda de sobrerrelaxação
        tolerancia = 0.0001; % criterio de parada do erro
        limite_iters = 100000; % numero maximo de iterações
        iter = 0; % contador de iteracoes
        erro_max = inf; % erro inicial para entrar no loop

        % Calcula A por liebmann até alcançar a tolerância desejada ou bater o
        % limite de iterações
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

                    % Escolhe qual equação de Aij utilizar
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
    
    elseif item == "e2"
        
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
                
        instantes = tempos(end)-1;
        Tempo = 0:dt:instantes*dt; % Vetor de tempos
        iter = 0; % contador de iteracoes
        
        Ainstantes = repmat(A, [1 1 3]);
        A = repmat(A, [1 1 2]); % Guarda o estado atual e anterior da matriz A
        k = 1;
        
        lambda = 1.75; % lambda de sobrerrelaxação
        
        for t=Tempo
            
            A(: , :, k) = A(:, :, k+1);
            iter = iter + 1;
            
            % Itera bobina por linha e coluna
            for j=rows_bobina 
                for i=cols_bobina
                    
                    % Se está embaixo ou em cima da bobina externa, pula linha
                    if Fronteiras(j, i) == vazio_direita
                        break;
                    end
                    
                    if any([row_eq(4) row_eq(16)] == j) && ...
                            any([col_eq(14) col_eq(16) col_eq(20)] == i)
                        continue
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
                    else
                        A(j, i, k+1) = Aij_int_bobina( ...
                            A, j, i, k, t, ...
                            MIx(j, i), MIy(j, i), JZ(j, i), Sigma(j, i));                                   
                    end
                        
                end
            end
             
            Ak = A(:, :, k+1);
            
            for j=2:rows-1
                for i=2:col_eq(20)-1

                    % Se for um ponto da bobina pula 
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

                    % Calcula valor novo de A no ponto por Liebmann
                    A(j, i, k+1) = lambda*Acalc + (1-lambda)*A(j, i, k);

                end
            end
           
            for idx=1:length(tempos)    
                if iter == tempos(idx)
                    Ainstantes(:, :, idx) = A(:, :, k+1);
                end
            end
            
        end
        
        A = Ainstantes;
        
    end
   
end

