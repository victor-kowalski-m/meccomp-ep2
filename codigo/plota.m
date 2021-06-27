function [] = plota(A, dx, dy, Bx, By, Hx, Hy, tempos)
    
    [X,Y] = meshgrid(0:dx:0.22,-0.1:dy:0.1);
    tempo = "";
    idx = 1;
    
    for k=1:length(tempos)
        
        if length(tempos) > 1
            tempo = " - tempo " + (tempos(k)-1) + " dt";
        end
        
        % Plot de Az
        figure("Name", "Az" + tempo);
        surf(X,Y,A(:, :, tempos(k)), 'LineStyle', ':')
        colorbar
        xlabel('x (m)')
        ylabel('y (m)')
        zlabel('Az')

%         % Plot de B
%         figure("Name", "B" + tempo);
%         quiver(X,Y,Bx(:, :, idx),By(:, :, idx));
%         axis equal
%         xlabel("x(m)")
%         ylabel("y(m)")
%         zlabel('B')
%         grid()
% 
%         % Plot de H
%         figure("Name", "H" + tempo);
%         quiver(X,Y,Hx(:, :, idx),Hy(:, :, idx));
%         axis equal
%         xlabel("x(m)")
%         ylabel("y(m)")
%         zlabel('H')
%         grid()
    
        idx = idx + 1;
        
    end
        
end

