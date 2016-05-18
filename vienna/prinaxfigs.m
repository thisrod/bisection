
# plot section through centre of mass
[~,i0] = min(abs(z1));  [X,Y] = ndgrid(x,y);
clf, hold on, axis([-4 4 100 106.5]), axis equal
contour(x, y, K(:,:,i0).', (1:9).^2/0.288.^2), colormap([0 0 0])
plot(X, Y, '+k', 'markersize', 2, 'color', [0.7 0.7 0.7])
plot(r0(1), r0(2), 'ok', 'markersize', 4)
text(0.05, 103.35, 'r_0', 'fontsize', 14)
plot(1.35+103.15i+0.288*exp(2*pi*1i*(0:0.01:1)), '-r')
text(1.4, 102.7, 'x_0', 'color', 'r', 'fontsize', 14)
xlabel x,  ylabel y
title(['Section of split trap, t = 17 ms, z = ' num2str(z(i0), 3)], 'fontsize', 14)
text(1.5, 105.9, 'K^2=1/x_0^2, 4/x_0^2, 9/x_0^2, ...', 'fontsize', 14)
