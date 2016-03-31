% load free gas modes

load debug.mat

gs = r.points(2);

norm(imag(ew)) / norm(ew), norm(imag(ev), 'fro') / norm(ev, 'fro'), norm(imag(U), 'fro') / norm(U, 'fro')

[~, i] = sort(abs(ew));  ew = real(ew(i));  ev = ev(:,i);
bmod = real(U*ev);  ev = real(ev);
buv = reshape(bmod,gs,2,[]);

figure, plot(r.x, LAP*a, r.x, K*a, r.x, r.a.g*diag(a.^2)*a)

figure, semilogy(abs(ew), '.k'), title Eigenvalues
% w = ck, k is around n/4L, healing length is (2*gamma)^{-1/2} in LL normalisation
% CHECKME
n = 1:length(ew);  L = r.ranges(2);  L = L*gs/(gs-1);
hold on, plot(n, (n-2)*sqrt(r.a.gamma)/(4*L), '-k')
legend computed sound Location SouthEast

figure, imagesc(bmod), colormap gray, title Eigenvectors

rsdl = sqrt(sum(abs(BdG*bmod-bmod*diag(ew)).^2));
ix = 1:numel(ew);  iy = ix(ew>0);  iz = ix(ew<0);
figure, semilogy(rsdl, '.k'), title 'Residuals of eigenvalue problem'

nms = squeeze(min(sqrt(sum(abs(reshape(buv, prod(gs), 2, []).^2)))));
figure, semilogy(nms, '.k'), title 'Eigenvector norms'

figure, plot(r.x, r.a.g*a.^2-mu, '-k', r.x, r.a.K(r), '-r')
legend('\mu discrep.', 'K')

k1 = 2*pi/L;

[sv,sw] = eig(-LAP, 'vector');
n = 1:length(sw);  k = k1*n/2;
figure, plot(n, sw, '.k', n, k.^2, '-k'), title('laplacian eigenvalues')
figure, imagesc(sv), colormap gray, title('laplacian eigenvectors')

[osv,osw] = eig(-U1'*LAP*U1, 'vector');
[~, i] = sort(abs(osw));  osw = osw(i);  osv = osv(:,i);
figure, plot(real(osw), '.k'), title('orthogonal laplacian eigenvalues')
figure, imagesc(real(U1*osv)), colormap gray, title('orthogonal laplacian eigenvectors')

figure
k1 = 2*pi/L;  k2 = 10*k1;
plot(r.x, sin(k1*r.x'), '-k', r.x, -k1^-2*LAP*sin(k1*r.x'), '--r')
hold on, 
plot(r.x, sin(k2*r.x'), '-k', r.x, -k2^-2*LAP*sin(k2*r.x'), '--r')
title 'sine waves reproduced with spectral Laplacian'

figure
plot(r.x, sv(:,2:3), '-k', r.x, -k1^-2*LAP*sv(:,2:3), '--r', r.x, -LAP*sv(:,2:3)/sw(2), ':b')
title 'numerical laplacian eigenvectors'

figure, plot(max(abs(LAP*sv(:,2:end)/diag(sw(2:end))+sv(:,2:end))), '.k')
title 'Laplacian eigenproblem residuals'