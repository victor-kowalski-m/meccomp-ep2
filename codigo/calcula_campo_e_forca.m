function [Bx, By, Hx, Hy, Fela_x, Fela_y]...
    = calcula_campo_e_forca(A, mi0, MIx, MIy, dx, dy, rows, cols, ...
    col_eq, Fronteiras, vazio_direita)
    
    % Inicialização das matrizes de densidade superficial de fluxo magnético e
    % intensidade de campo magnético
    [Bx, By, Hx, Hy] = deal(zeros(rows, cols));

    for j = 2:rows-1
       for i = 2:cols-1

            % Se está embaixo ou em cima da bobina externa, pula linha
            if Fronteiras(j, i) == vazio_direita
                break;
            end

           Bx(j,i) = (A(j+1,i) - A(j-1,i))/(2*dy);
           By(j,i) = -(A(j,i+1) - A(j,i-1))/(2*dx);
           Hx(j,i) = Bx(j,i)/MIx(j,i);
           Hy(j,i) = By(j,i)/MIy(j,i);

       end
    end
    
    Fela_x = 1/(2*mi0)*(...
    trapz(-0.1:dx:0.1, Bx(:, col_eq(4)).^2 - By(:, col_eq(4)).^2) );

    Fela_y = 1/(2*mi0)*(...
        trapz(-0.1:dx:0.1, 2*Bx(:, col_eq(4)).*By(:, col_eq(4))) );


    fprintf('A força Fx = %8.4f \n',Fela_x);
    fprintf('A força Fy = %8.4f \n',Fela_y);

end

