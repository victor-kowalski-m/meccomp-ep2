function [figA, figB, figH] = plota(A, dx, dy, Bx, By, Hx, Hy)
    
    [X,Y] = meshgrid(0:dx:0.22,0:dy:0.20);
    
    % Plot de Az
    figA = figure("Name", "Vetor Az");
    surf(X,Y,A, 'LineStyle', ':')
    colorbar
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('Az')

    % Plot de B
    figB = figure("Name", "Vetor de densidade superficial de fluxo magnético");
    quiver(X,Y,Bx,By);
    axis equal
    xlabel("x(m)")
    ylabel("y(m)")
    zlabel('B')
    grid()

    % Plot de H
    figH = figure("Name", "Vetor de intensidade de campo magnético");
    quiver(X,Y,Hx,Hy);
    axis equal
    xlabel("x(m)")
    ylabel("y(m)")
    zlabel('H')
    grid()

end

