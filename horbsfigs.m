set(0, 'defaultaxesfontsize', 14, 'defaulttextfontsize', 14)
load horbsplotdat.mat

l = [length(x), length(y)];

for t = [0 12:17], for j = [2 size(w,2)], for n = 0:1

thisw = w(t+1,j,:,n+1);
figure, zplot(x,y, reshape(thisw,l)'), axis equal, axis off, hold on
contour(x, y, reshape(K(t+1,:), l)', (1:9).^2/0.25.^2)
colormap(0.7*[1 1 1])
if t == 0 && n == 0
	patch([-2.7 -2.7 -1.7 -1.7], [-0.4 -0.5 -0.5 -0.4], 0.5*[1 1 1])
	text(-2.5, -0.3, '1 \mum', 'color', 'w')
end
if j == 2
	s = sprintf('qmode%d-%d.pdf', t, n);
else
	s = sprintf('wmode%d-%d.pdf', t, n);
end

saveTightFigure(s)

end, end, end