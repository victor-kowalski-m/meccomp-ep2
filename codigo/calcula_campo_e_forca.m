function [Bx, By, Hx, Hy, Fela_x, Fela_y]...
    = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, cols, ...
    col_eq, Fronteiras, vazio_direita, k, printa)
    
    % Inicializa matrizes vazias de B e H
    [Bx, By, Hx, Hy] = deal(zeros(rows, cols));
    
    % Calcula B e H em cada ponto do domínio
    for j = 2:rows-1
       for i = 2:cols-1

            % Se está embaixo ou em cima da bobina externa, pula linha
            if Fronteiras(j, i) == vazio_direita
                break;
            end 
           
            % Primeira diferença central para B
            Bx(j,i) = (A(j+1,i,k) - A(j-1,i,k))/(2*dy);
            By(j,i) = -(A(j,i+1,k) - A(j,i-1,k))/(2*dx);
            
            % H é a divisão de B no ponto pelo mi no ponto
            Hx(j,i) = Bx(j,i)/MIx(j,i);
            Hy(j,i) = By(j,i)/MIy(j,i);

       end
    end

    % Cálcula Fela em x e y por integral de linha (método dos trapézios)
    Fela_x = 1/(2*mi0)*(...
    trapz(-0.1:dx:0.1, Bx(:, col_eq(4)).^2 - By(:, col_eq(4)).^2) );
    Fela_y = 1/(2*mi0)*(...
        trapz(-0.1:dx:0.1, 2*Bx(:, col_eq(4)).*By(:, col_eq(4))) );

    % Se o argumento "printa" foi setado, imprime valores de Fela
    if printa
        fprintf("Fela_x = %8.4f \n",Fela_x);
        fprintf("Fela_y = %8.4f \n",Fela_y);
    end
    
end

