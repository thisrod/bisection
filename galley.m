[X,Y,Z] = ndgrid(x,y,z);  r = [X(:) Y(:) Z(:)]-repmat(r0', numel(X), 1);
r = r.';  Kfit = reshape(sum(r.*(U*r)), size(K));

[~,j] = min(abs(z-r0(3)));   WHITE = [1 1 1];
figure, hold on
plot(X(:), Y(:), '+', r0(1), r0(2), 'o', 'color', 0.65*WHITE)
contour(x, y, K(:,:,j).', (0:0.2:2).^2, '-r')
% contour(x, y, K(:,:,j)-Kfit(:,:,j), (0:0.2:2).^4, '-k')