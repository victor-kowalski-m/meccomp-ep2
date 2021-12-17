function [] = plota(A, dx, dy, Bx, By, Hx, Hy, tempos, Fela_x, Fela_y, dt)
    
    [X,Y] = meshgrid(0:dx:0.22,-0.1:dy:0.1); % matrizes com valores de X e Y em todo o domínio
    tempo = ""; % string auxiliar para indicar instante de tempo do plot
    
    % Itera cada instante de tempo calculado
    for k=1:length(tempos)
        
        % Se houver mais de um tempo, indica no título do plot
        if length(tempos) > 1
            tempo = " " + (tempos(k)) + "dt";
        end
        
        % Plot de Az
        figure("Name", "Az" + tempo, 'NumberTitle','off');
        surf(X,Y,A(:, :, k), 'LineStyle', ':')
        colorbar
        xlabel('x (m)')
        ylabel('y (m)')
        zlabel('Az')

        % Plot de B
        figure("Name", "B" + tempo, 'NumberTitle','off');
        quiver(X,Y,Bx(:, :, k),By(:, :, k));
        axis equal
        xlabel("x(m)")
        ylabel("y(m)")
        zlabel('B')
        grid()

        % Plot de H
        figure("Name", "H" + tempo, 'NumberTitle','off');
        quiver(X,Y,Hx(:, :, k),Hy(:, :, k));
        axis equal
        xlabel("x(m)")
        ylabel("y(m)")
        zlabel('H')
        grid()
        
    end
    
    % Se houver mais de um instante de tempo calculado, plota Fela
    if length(tempos) > 1
        Tempo = 0:dt:tempos(end)*dt; % eixo horizontal de tempo
        
        % Plot de Fela_x
        figure("Name", "Fela_x", 'NumberTitle','off');
        plot(Tempo, Fela_x);
        grid();
        xlim([Tempo(1) Tempo(end)]);
        
        % Plot de Fela_x
        figure("Name", "Fela_y", 'NumberTitle','off');
        plot(Tempo, Fela_y);
        grid();
        xlim([Tempo(1) Tempo(end)]);
    end
        
end

