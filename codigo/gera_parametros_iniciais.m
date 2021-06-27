function ...
    [dx, dy, rows, cols, col_eq, row_eq, A, JZ, Fronteiras, vertical, ...
     horizontal, vazio_direita, mi0, MIx, MIy, Sigma, dt, tempos]...
     = gera_parametros_iniciais(item, dx, dt)
    
    if ~any(item == ["ad" "e1" "e2"])
        disp("Item escolhido é inválido!")
        return
    end
 
    [dx, dy] = deal(dx); % passos em metros

    rows = 0.20/dy + 1; % numero de linhas da matriz
    cols = 0.22/dx + 1; % numero de colunas da matriz

    % Funções para calcular os pontos da malha equivalentes a x e y no desenho
    col_eq = @(x) round(x/22*(cols-1)+1); 
    row_eq = @(y) round(y/20*(rows-1)+1);
    
    A = zeros(rows, cols); % matriz inicial com zeros
    
    % Define matriz de e densidade superficial de corrente elétrica no domínio
    JZ = zeros(rows, cols);
    Jz = @(y) (2*10^6*cos(pi*(y-0.1)/(12*10^-2)) + 8*10^5);
    for j=row_eq(4)+1:row_eq(16)-1 
        for i = col_eq(20)+1:col_eq(22)-1
            JZ(j, i) = Jz((j-1)*dy);
        end
        for i = col_eq(14)+1:col_eq(16)-1 
            JZ(j, i) = -Jz((j-1)*dy);  
        end
    end

    % Define fronteiras do domínio: 0 interior, 1 vertical, 2 horizontal
    vertical = 1;
    horizontal = 2;
    vazio_direita = 3;
    Fronteiras = zeros(rows, cols);
    Fronteiras(2:end-1, col_eq(4)) = vertical;
    Fronteiras([2:row_eq(4) row_eq(16):end-1], col_eq(5)) = vertical;
    Fronteiras(row_eq(4):row_eq(16), [col_eq(16) col_eq(20)]) = vertical;
%     Fronteiras(row_eq(4)+1:row_eq(16)-1, col_eq(14)) = vertical;
    Fronteiras([row_eq(4) row_eq(16)], col_eq(5):col_eq(16)) = horizontal;
    Fronteiras([2:row_eq(4) row_eq(16):end-1], col_eq(20):col_eq(22)-1)...
        = vazio_direita;
    
    mi0 = 4*pi*10^-7; % permeabilidade magnética do vacuo, do ar e da bobina
    Sigma = 0;
    tempos = 1;
    
    if item == "ad"
        miferro = 2500*mi0; % permeabilidade magnética do ferro
        MI = ones(rows, cols)*mi0;
        MI(:, col_eq(0):col_eq(4)) = miferro;
        MI(:, col_eq(16):col_eq(20)) = miferro;
        MI(row_eq(0):row_eq(4), col_eq(5):col_eq(20)) = miferro;
        MI(row_eq(16):row_eq(20), col_eq(5):col_eq(20)) = miferro;
        
        MIx = MI;
        MIy = MI;
          
    else
        % Define matriz de permeabilidade magnética no domínio
        mi0 = 4*pi*10^-7; % permeabilidade magnética do vacuo, do ar e da bobina
        miferrox = 1200*mi0; % 1200 permeabilidade magnética do ferro em x
        miferroy = 2500*mi0; % permeabilidade magnética do ferro em y
        [MIx,MIy] = deal(ones(rows, cols)*mi0);
        MIx(:, col_eq(0):col_eq(4)) = miferrox;
        MIx(:, col_eq(16):col_eq(20)) = miferrox;
        MIx(row_eq(0):row_eq(4), col_eq(5):col_eq(20)) = miferrox;
        MIx(row_eq(16):row_eq(20), col_eq(5):col_eq(20)) = miferrox;
        MIy(:, col_eq(0):col_eq(4)) = miferroy;
        MIy(:, col_eq(16):col_eq(20)) = miferroy;
        MIy(row_eq(0):row_eq(4), col_eq(5):col_eq(20)) = miferroy;
        MIy(row_eq(16):row_eq(20), col_eq(5):col_eq(20)) = miferroy;
        
        if item == "e1"
            
        elseif item == "e2"
            sigma = 4*10^6;
            Sigma = zeros(rows, cols);
            Sigma(row_eq(4):row_eq(16),...
                [col_eq(14):col_eq(16) col_eq(20):col_eq(22)])...
                = sigma;
            
            tempos = [10 100 500];
        end 
    end
end

