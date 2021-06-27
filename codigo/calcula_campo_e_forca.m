function [Bx, By, Hx, Hy, Fela_x, Fela_y]...
    = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, cols, ...
    col_eq, Fronteiras, vazio_direita, k, printa)
    
    % Inicialização das matrizes de densidade superficial de fluxo magnético e
    % intensidade de campo magnético
    [Bx, By, Hx, Hy] = deal(zeros(rows, cols));
        
    for j = 2:rows-1
       for i = 2:cols-1

            % Se está embaixo ou em cima da bobina externa, pula linha
            if Fronteiras(j, i) == vazio_direita
                break;
            end 

           Bx(j,i) = (A(j+1,i,k) - A(j-1,i,k))/(2*dy);
           By(j,i) = -(A(j,i+1,k) - A(j,i-1,k))/(2*dx);
           Hx(j,i) = Bx(j,i)/MIx(j,i);
           Hy(j,i) = By(j,i)/MIy(j,i);

       end
    end

    Fela_x = 1/(2*mi0)*(...
    trapz(-0.1:dx:0.1, Bx(:, col_eq(4)).^2 - By(:, col_eq(4)).^2) );

    Fela_y = 1/(2*mi0)*(...
        trapz(-0.1:dx:0.1, 2*Bx(:, col_eq(4)).*By(:, col_eq(4))) );

    if printa
        fprintf("Fela_x = %8.4f \n",Fela_x);
        fprintf("Fela_y = %8.4f \n",Fela_y);
    end
    
end

