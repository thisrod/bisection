set(0, 'defaultaxesfontsize', 14, 'defaulttextfontsize', 14)
load horbsplotdat.mat

for t = [0 12:17], for n = 0:1

l = [length(x), length(y)];
figure, zplot(x,y, reshape(q(t+1,:,n+1), l).'), axis equal, axis off, hold on
contour(x, y, reshape(K(t+1,:), l)', (1:9).^2/0.25.^2)
colormap(0.7*[1 1 1])
if t == 0 && n == 0
	patch([-2.7 -2.7 -1.7 -1.7], [-0.4 -0.5 -0.5 -0.4], 0.5*[1 1 1])
	text(-2.5, -0.3, '1 \mum', 'color', 'w')
end
s = sprintf('qmode%d-%d.pdf', t, n);
saveTightFigure(s)

end, end