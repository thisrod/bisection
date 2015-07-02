load horbsplotdat.mat
clf
zplot(x,y,ql'), axis equal, axis off, hold on
contour(x, y, reshape(K, length(x), length(y))', (1:9).^2/0.288.^2)
colormap(0.7*[1 1 1])
print('-depsc', ['qmode' int2str(t) '-0'])