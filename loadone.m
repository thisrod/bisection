% load free gas modes

load debug1.mat

assert(isa(U, 'double'))

fprintf('DEBUGGING INFORMATION FOR 1D BdG MODES\n\n')
fprintf('Machine epsilon: %.2e\n\n', eps)
fprintf('Imaginary fraction of eigen-\n     values: %.2e\n     vectors: %.2e\n\n', ...
	norm(imag(ew)) / norm(ew), ...
	norm(imag(ev), 'fro') / norm(ev, 'fro'))

ixs = abs(real(ew)) > abs(imag(ew));
ewr = ew(ixs);  [~,j] = sort(abs(ewr));  ewr = ewr(j);
ewi = ew(~ixs);  [~,j] = sort(abs(ewi), 'descend');  ewi = ewi(j);
figure, plot(real([ewi; ewr]), '.k'), hold on, plot(imag([ewi; ewr]), '.r')
title(sprintf('Re and Im of ew, mu = %.8f', mu))

assert(norm(imag(ew)) / norm(ew) < 1e-5), ew = real(ew);
pew = ew(ew>0);  new = ew(ew<0);
fprintf('Mismatch between +ive and -ive eigenvalues: %.2e\n\n', ...
	norm(sort(pew) + sort(new, 'descend')))
	
fprintf('non-orthonormality of eigenvectors: %.2e\n', ...
	norm(ev'*ev-eye(2*gs-2), 'fro'))

% remove u/v degeneracy and sort by increasing eigenvalue
% note that sin/cos degeneracy remains
pev = ev(:, ew>0);  [~, i] = sort(pew);
ew = pew(i);  ev = pev(:,i);

evIm = sqrt(sum(imag(ev).^2));
figure, plot(evIm, '.k')
title('Imaginary components of BdG eigenvectors')
ylabel('|Im u_i|^2')

% Fix Matlab quirk where four modes are elliptically polarised
bmod = U*ev;
ixodd = find(evIm > 0.1);
for i = 1:2:length(ixodd)
	ix = ixodd(i);
	bmod(:,ix:ix+1) = [real(bmod(:,ix)) imag(bmod(:,ix))];
	bmod(:,ix) = bmod(:,ix) / norm(bmod(:,ix));
	bmod(:,ix+1) = bmod(:,ix+1) / norm(bmod(:,ix+1));
end
buv = reshape(bmod,gs,2,[]);

load debug2.mat

fprintf('non-orthonormality of corrected eigenvectors: %.2e\n', ...
	norm(bmod'*bmod-eye(gs-1), 'fro'))

figure
nodd = nnz(ixodd);  odd = buv(:, :, ixodd);  wodd = ew(ixodd);
for i = 1:nodd, for j = 1:2
	subplot(nodd, 2, 2*i+j-2)
	plot(r.x, real(odd(:,j,i)), '-k', r.x, imag(odd(:,j,i)), '-r')
	title(num2str(wodd(i)))
end, end

% Eigenstuff is real below here
bmod = real(bmod);  buv = real(buv);  ev = real(ev);

% figure, plot(r.x, LAP*a, r.x, K*a, r.x, r.a.g*diag(a.^2)*a)

figure, plot(ew, '.k'), title('computed and expected BdG eigenvalues')
% w = ck, k is around n/2L, healing length is (2*gamma)^{-1/2} in LL normalisation
% FIXME include ordinary energy
n = 1:length(ew);  L = r.ranges(2);  L = L*gs/(gs-1);
k1 = 2*pi/L;  k = k1*n/2;  heal = 1/sqrt(2*r.a.gamma);
hold on, plot(n, k.*sqrt(heal.^-2+k.^2), '-k')
% plot(n, (n-2)*sqrt(r.a.gamma)/(4*L), '-k')

figure, subplot 211, imagesc(squeeze(buv(:,1,:))), title 'BdG u modes'
subplot 212, imagesc(squeeze(buv(:,2,:))), title 'BdG v modes'
xlabel mode, ylabel x
colormap gray

rsdl = sqrt(sum(abs(BdG*bmod-bmod*diag(ew)).^2));
ix = 1:numel(ew);  iy = ix(ew>0);  iz = ix(ew<0);
figure, semilogy(rsdl, '.k'), title 'Residuals of eigenvalue problem'

nms = squeeze(sqrt(2*sum(abs(buv).^2)));
nms = nms - 1;
figure, plot(n, nms(1,:), '.k', n, nms(2,:), '.r'), legend u v
title '1D Bogoliubov coefficients minus 1'

u1 = squeeze(buv(:,1,1));  v1 = squeeze(buv(:,2,1));
figure, subplot 311, plot(r.x, u1, '-k', r.x, v1, '-r')
title('Lowest BdG u and v modes')
subplot 312, plot(r.x, -LAP*u1, r.x, r.a.K(r)'.*u1, r.x, -mu*u1, ...
	r.x, 2*r.a.g*a.^2.*u1, r.x, -r.a.g*a.^2.*v1)
legend lap trap mu self other 
subplot 313, plot(r.x, (-Bself+M)*u1+Bother*v1, ':k', r.x, ew(1)*u1, '--k')
title 'BdG equation CHECKME', legend LHS '\epsilon u_1'

figure, plot(r.x, r.a.g*a.^2-mu, '-k', r.x, r.a.K(r), '-r')
legend('\mu discrep.', 'K')

return


[sv,sw] = eig(-LAP, 'vector');
n = 1:length(sw);  k = k1*n/2;
figure, plot(n, sw, '.k', n, k.^2, '-k'), title('laplacian eigenvalues')
figure, imagesc(sv), colormap gray, title('laplacian eigenvectors')

[osv,osw] = eig(-U1'*LAP*U1, 'vector');
[~, i] = sort(abs(osw));  osw = osw(i);  osv = osv(:,i);
figure, plot(real(osw), '.k'), title('orthogonal laplacian eigenvalues')
figure, imagesc(real(U1*osv)), colormap gray, title('orthogonal laplacian eigenvectors')

figure
k2 = 10*k1;
plot(r.x, sin(k1*r.x'), '-k', r.x, -k1^-2*LAP*sin(k1*r.x'), '--r')
hold on, 
plot(r.x, sin(k2*r.x'), '-k', r.x, -k2^-2*LAP*sin(k2*r.x'), '--r')
title 'sine waves reproduced with spectral Laplacian'

figure
plot(r.x, sv(:,2:3), '-k', r.x, -k1^-2*LAP*sv(:,2:3), '--r', r.x, -LAP*sv(:,2:3)/sw(2), ':b')
title 'numerical laplacian eigenvectors'

figure, plot(max(abs(LAP*sv(:,2:end)/diag(sw(2:end))+sv(:,2:end))), '.k')
title 'Laplacian eigenproblem residuals'