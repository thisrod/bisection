% load BdG modes and process

load debug.mat

a = a0(:);

gs = r.points(2:end);
bop = [-LAP+NL, -r.c.rpsn*diag(a.^2);
	r.c.rpsn*diag(conj(a).^2), LAP-NL];

assert(norm(imag(ew)) < 1e-10*norm(ew))
assert(norm(imag(ev), 'fro') < 1e-2*norm(ev, 'fro'))
assert(norm(imag(U), 'fro') < 1e-2*norm(U, 'fro'))

[~, i] = sort(abs(ew));  ew = real(ew(i));  ev = ev(:,i);
bmod = real(U*ev);  ev = real(ev);
buv = reshape(bmod,gs(1),gs(2),2,[]);

figure, semilogy(abs(ew), '.k')

figure
for j = 1:2, subplot(2,1,j), surf(buv(:,:,j,1)), end

rsdl = sqrt(sum(abs(bop*bmod-bmod*diag(ew)).^2));
ix = 1:numel(ew);  iy = ix(ew>0);  iz = ix(ew<0);
figure, semilogy(rsdl, '.k');

nms = squeeze(min(sqrt(sum(abs(reshape(buv, prod(gs), 2, []).^2)))));
figure, semilogy(nms, '.k')