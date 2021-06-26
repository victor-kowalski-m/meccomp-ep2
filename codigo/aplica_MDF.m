function [A, iters, erro_max] = ...
    aplica_MDF(item, A, rows, cols, Fronteiras, vertical, horizontal, ...
    vazio_direita, Aij_interior, Aij_front_vert, Aij_front_hori, ...
    MIx, MIy, JZ, lenT)
    
    if any(item == ["ad" "e1"])
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
                        Acalc = Aij_front_vert(A, j, i, ...
                        MIx(j,i-1), MIy(j,i-1), JZ(j,i-1), ...
                        MIx(j,i+1), MIy(j,i+1), JZ(j,i+1));
                    elseif Fronteiras(j, i) == horizontal % se é fronteira horizontal
                        Acalc = Aij_front_hori(A, j, i, ...
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
        
        T = 0:dt:(lenT-1)*dt; % Vetor de tempos
        A = repmat(A0,[1 1 lenT]); % Matriz inicial de valores de A no tempo

        rows_bobina = row_eq(4):row_eq(16);
        cols_bobina = [col_eq(14):col_eq(16) col_eq(20):col_eq(22)];

        for k = 1:lenT - 1

            % Itera bobina por linha e coluna
            for j=rows_bobina 
                for i=cols_bobina

                    % Escolhe qual equação de Aij utilizar
                    if MIx(j, i-1) ~= MIx(j, i+1) % se é fronteira vertical
                        A(j, i, k+1) = Aijkmais1_front_vert( ...
                            A, j, i, k, T(k), ...
                            MIx(j, i-1), MIy(j, i-1), JZ(j,i-1), ...
                            MIx(j, i+1), MIy(j, i+1), JZ(j,i+1), sigma1, sigma2);

                    elseif MIx(j-1,i) ~= MIx(j+1,i) % se é fronteira horizontal
                        A(j, i, k+1) = Aijkmais1_front_hori( ...
                            A, j, i, k, T(k), ...
                            MIx(j, i-1), MIy(j, i-1), JZ(j,i-1), ...
                            MIx(j, i+1), MIy(j, i+1), JZ(j,i+1), sigma1, sigma2);

                    else % se estiver em um ponto interior do domínio
                        A(j, i, k+1) = Aijkmais1_interior( ...
                            A, j, i, k, T(k), ...
                            MIx(j, i), MIy(j, i), JZ(j, i), sigma1, sigma2);

                    end                                    

                end
            end

            % Itera pontos fora da bobina

        end
    end
   
end

