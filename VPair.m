% Find a bogoliubov ground state for the Vienna trap

op = static(vienna, 0);
op.ranges(1) = 0.5;
op.points(1) = 49;
op.steps = 30;
op = eqop(op);

uvs = bogs(op);
U = real(uvs.a.U);  V = real(uvs.a.V);

figure, plot(uvs.a.bew, '.k')
for i = [1 2 5 25 100 200 500 1e3]
	figure, subplot 211, imagesc(reshape(U(:,i), uvs.points(2:end)))
	subplot 212, imagesc(reshape(V(:,i), uvs.points(2:end))), colormap gray
end