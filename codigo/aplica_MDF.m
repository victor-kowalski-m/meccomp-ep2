function [A, iters] = ...
    aplica_MDF(item, A, rows, cols, row_eq, col_eq, Fronteiras, ...
    vertical, horizontal, vazio_direita, equacoes, MIx, MIy, JZ, Sigma, dt)
    
    if any(item == ["ad" "e1"])
        
        Aij_interior = equacoes{1};
        Aij_vert = equacoes{2};
        Aij_hori = equacoes{3};
        
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
        
%         % Equações de pontos interiores fora da bobina
%         Aij_int_esquerda = equacoes{4};
%         Aij_int_direita = equacoes{5};
%         Aij_int_abaixo = equacoes{6};
%         Aij_int_acima = equacoes{7};
%         
%         % Equações de fronteira vertical fora da bobina
%         Aij_vert_esquerda = equacoes{8};
%         Aij_vert_direita = equacoes{9};
%         Aij_vert_abaixo = equacoes{10};
%         Aij_vert_acima = equacoes{11};
%         
%         % Equações de fronteira horizontal fora da bobina
%         Aij_hori_esquerda = equacoes{12};
%         Aij_hori_direita = equacoes{13};
%         Aij_hori_abaixo = equacoes{14};
%         Aij_hori_acima = equacoes{15};
        
        instantes = 500;
        Tempo = 0:dt:instantes*dt; % Vetor de tempos
        A = repmat(A,[1 1 instantes+1]); % Matriz inicial de valores de A no tempo
        rows_bobina = row_eq(4):row_eq(16);
        cols_bobina = [col_eq(14):col_eq(16) col_eq(20):col_eq(22)-1];
        iters = 0; % contador de iteracoes
        
        for k = 1:instantes
            
            iters = iters + 1;
            
            % Itera bobina por linha e coluna
            for j=rows_bobina 
                for i=cols_bobina

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
                            MIx(j, i-1), MIy(j, i-1), JZ(j,i-1), ...
                            MIx(j, i+1), MIy(j, i+1), JZ(j,i+1), ...
                            Sigma(j, i-1), Sigma(j, i+1));

                    else % se estiver em um ponto interior do domínio
                        A(j, i, k+1) = Aij_int_bobina( ...
                            A, j, i, k, Tempo(k), ...
                            MIx(j, i), MIy(j, i), JZ(j, i), Sigma(j, i));
                    end                                    

                end
            end
            
            % Itera pontos à esquerda da(s) bobinas
            

        end
    end
   
end

