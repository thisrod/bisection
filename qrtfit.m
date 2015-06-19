% fit quadratic/quartic potentials to the Vienna trap data, for interpolation around the bottom of the trap.

set(0, 'defaultaxesfontsize', 14, 'defaulttextfontsize', 14)

t = 13;  prinax,  [x, y, z, K] = loadk(t);
[r0, wgt] = tcent(x,y,z,K);
x1 = x - r0(1);  y1 = y - r0(2);  z1 = z - r0(3);

% fit symmetric quartic to wells below second excited state
[~,i0] = min(abs(y1));  [~,j0] = min(abs(z1));
K1 = squeeze(K(:,i0,j0));
[~,i1] = min(abs(x1));  ft = (K1 <= 50);
xf = x(ft)';  kf = K1(ft);
xcfs = [xf.^4 xf.^2 ones(size(xf))] \ kf;
K2 = [x'.^4 x'.^2 ones(size(x'))]*xcfs;

% r expands sample points over principal axes
[X,Y,Z] = ndgrid(x1,y1,z1);  r = [X(:) Y(:) Z(:)]*U;
rsdl = K(:) - [r(:,1).^4 r(:,1).^2 ones(numel(X),1)]*xcfs;
cfs = ([r(:,2).^2 r(:,3).^2].*repmat(wgt,1,2)) \ (rsdl.*wgt);
Kfit = [r(:,1).^4 r(:,1).^2 ones(numel(X),1) r(:,2).^2 r(:,3).^2]*[xcfs; cfs];
Kfit = reshape(Kfit, size(K));

figure, plot(x,K1,'-k', x, K2, '-r')
title 'Fit along splitting axis'
text(-1.5, 125, sprintf('%.2f x^4 %+.2f x^2 %+.2f', xcfs(1), xcfs(2), xcfs(3)), 'color', 'r')
xlabel x, ylabel K^2
axis([-3 3 0 200])

figure
for i = 1:2
	r1 = r0+U*[0 0 z(j0)]';
	subplot(1,2,i), hold on
	contour(x, y, K(:,:,j0).', (1:10).^2/0.288.^2), colormap([0 0 0])
	contour(x, y, Kfit(:,:,j0).', (1:10).^2/0.288.^2, '-r')
	plot(r1(1), r1(2), 'ok', 'markersize', 4)
	xlabel x, ylabel y, axis equal, axis([-4 4 100 106])
	title(sprintf('Fit at t = %d, z = %.1f', t, z(j0)))
	j0 = 70;
endfor
