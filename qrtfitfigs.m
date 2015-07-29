set(0, 'defaultaxesfontsize', 14, 'defaulttextfontsize', 14)
load qrtfitdat

xx = linspace(-2.2, 2.2, 100)';
ctrs = (1:10).^2/0.25.^2;

for i = 1:length(tplot)

t = tplot(i);

figure, plot(x1(1,:), Ks1(i,:), '.k')
axis([-2.2 2.2 0 200]), hold on
plot(xx, [xx.^4 ones(size(xx)) xx.^2]*cfs(1:3, t+1), '-k')
title(['Fit along splitting axis at t = ' int2str(t) ' ms'])
s = sprintf('K^2 = %.2f x^4 %+.2f x^2 %+.2f', cfs(1, t+1), cfs(3, t+1), cfs(2, t+1));
text(-1.5, 125, s)
xlabel x, ylabel K^2

% TODO get Matlab to put a box around the figure

for j = 1:2
	figure
	contour(x1(j,:), y1, squeeze(Ks(i,:,:,j))', ctrs, '-k')
	hold on
	contour(x2, y2, squeeze(Kf(i,:,:,j))', ctrs, ':k')
	% xlabel x, ylabel y
	axis equal, axis([-4 4 -3 3]), axis off
	% title(sprintf('Fit at t = %d, z = %.1f', t, zsec(i,j)))
	if i == 1
		patch([-3.8 -3.8 -2.8 -2.8], [-1 -1.15 -1.16 -1], 0.5*[1 1 1])
		text(-3.5,-0.8,'1 \mum')
	end
	if j == 1, saveTightFigure(['kmode' int2str(t) '.pdf']), end
end

end	% i = 1:length(tplot)

figure
plot(0:17, cfs(1:4,:)), hold on
line(0:17, 1e4*cfs(5,:))
title 'coefficients of quartic fit', xlabel t
legend('x^4', '1', 'x^2', 'y^2', '10^4*z^2')
