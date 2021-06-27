function [] = plota(A, dx, dy, Bx, By, Hx, Hy, tempos, Fela_x, Fela_y, dt)
    
    [X,Y] = meshgrid(0:dx:0.22,-0.1:dy:0.1);
    tempo = "";
    idx = 1;
    
    for k=1:length(tempos)
        
        if length(tempos) > 1
            tempo = " " + (tempos(k)) + "dt";
        end
        
        % Plot de Az
        figure("Name", "Az" + tempo, 'NumberTitle','off');
        surf(X,Y,A(:, :, idx), 'LineStyle', ':')
        colorbar
        xlabel('x (m)')
        ylabel('y (m)')
        zlabel('Az')

        % Plot de B
        figure("Name", "B" + tempo, 'NumberTitle','off');
        quiver(X,Y,Bx(:, :, idx),By(:, :, idx));
        axis equal
        xlabel("x(m)")
        ylabel("y(m)")
        zlabel('B')
        grid()

        % Plot de H
        figure("Name", "H" + tempo, 'NumberTitle','off');
        quiver(X,Y,Hx(:, :, idx),Hy(:, :, idx));
        axis equal
        xlabel("x(m)")
        ylabel("y(m)")
        zlabel('H')
        grid()
    
        idx = idx + 1;
        
    end
    
    if length(tempos) > 1
        Tempo = 0:dt:tempos(end)*dt;
        figure("Name", "Fela_x", 'NumberTitle','off');
        plot(Tempo, Fela_x);
        grid();
        xlim([Tempo(1) Tempo(end)]);
        figure("Name", "Fela_y", 'NumberTitle','off');
        plot(Tempo, Fela_y);
        grid();
        xlim([Tempo(1) Tempo(end)]);
    end
        
end

