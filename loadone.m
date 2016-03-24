% load free gas modes

load debug.mat

gs = r.points(2);

norm(imag(ew)) / norm(ew), norm(imag(ev), 'fro') / norm(ev, 'fro'), norm(imag(U), 'fro') / norm(U, 'fro')

[~, i] = sort(abs(ew));  ew = real(ew(i));  ev = ev(:,i);
bmod = real(U*ev);  ev = real(ev);
buv = reshape(bmod,gs,2,[]);

figure, plot(r.x, LAP*a, r.x, K*a, r.x, r.a.g*diag(a.^2)*a)

figure, semilogy(abs(ew), '.k'), title Eigenvalues

figure, imagesc(bmod), colormap gray, title Eigenvectors

rsdl = sqrt(sum(abs(BdG*bmod-bmod*diag(ew)).^2));
ix = 1:numel(ew);  iy = ix(ew>0);  iz = ix(ew<0);
figure, semilogy(rsdl, '.k'), title 'Residuals of eigenvalue problem'

nms = squeeze(min(sqrt(sum(abs(reshape(buv, prod(gs), 2, []).^2)))));
figure, semilogy(nms, '.k'), title 'Eigenvector norms'

figure, plot(r.x, r.a.g*a.^2-mu, '-k', r.x, r.a.K(r), '-r')
legend('\mu discrep.', 'K')

[sv,sw] = eig(-LAP, 'vector');
figure, plot(sw, '.k'), title('laplacian eigenvalues')
figure, imagesc(sv), colormap gray, title('laplacian eigenvectors')

[osv,osw] = eig(-U1'*LAP*U1, 'vector');
[~, i] = sort(abs(osw));  osw = osw(i);  osv = osv(:,i);
figure, plot(real(osw), '.k'), title('orthogonal laplacian eigenvalues')
figure, imagesc(real(U1*osv)), colormap gray, title('orthogonal laplacian eigenvectors')