clf
for i = 1:2
	r1 = r0+U*[0 0 z(j0)]';
	subplot(1,2,i), hold on
	contour(x, y, K(:,:,j0).', (1:10).^2/0.288.^2), colormap([0 0 0])
	contour(x, y, Kfit(:,:,j0).', (1:10).^2/0.288.^2, '-r')
	plot(r1(1), r1(2), 'ok', 'markersize', 4)
	xlabel x, ylabel y, axis equal, axis([-4 4 100 106])
	title(sprintf('Quartic/quadratic fit, t = %d, z = %.1f', t, z(j0)))
	j0 = 70;
endfor
