% working script for Bogoliubov mode orthogonalisation

N = r.nspace -1;
modes = reshape(B,r.nspace,2,[]);
U = squeeze(modes(:,1,:));  V = squeeze(modes(:,2,:));
keys = round(ew, 10);		% identify ews that are within 1e-10

X = r.dV*(U'*U - V'*V);
nmsq = squeeze(r.dV*sum(abs(modes).^2));  nms = sqrt(nmsq);

new_modes = nan(size(modes));

for key = unique(keys)'
	tmp = find(keys == key);
	if numel(tmp) == 1
		new_modes(:,:,tmp) = modes(:,:,tmp);
		continue
	end
	assert(numel(tmp) == 2);
	j = tmp(1);  k = tmp(2);
	numerator = -X(j,k) + sqrt(X(j,k)^2-(nmsq(1,j)-nmsq(2,j))*(nmsq(1,k)-nmsq(2,k)));
	l = numerator/(nmsq(1,k)-nmsq(2,k));
	m = numerator/(nmsq(1,j)-nmsq(2,j));

	new_modes(:,:,j) = modes(:,:,j) + l*modes(:,:,k);
	new_modes(:,:,k) = modes(:,:,k) + m*modes(:,:,j);
end

c = r.dV*sum(abs(new_modes(:,1,:)).^2 - abs(new_modes(:,2,:)).^2);
c = 1./sqrt(c);
new_modes = new_modes.*repmat(c,r.nspace,2,1);

UU = squeeze(new_modes(:,1,:));  VV = squeeze(new_modes(:,2,:));

figure, plot(1:N, nms(1,:), 'ok', 1:N, nms(2,:), 'xk', [1:N; 1:N], nms, '-k' );
	title 'norms of U and V modes'
figure, plot(1:N, real(ew), '.k', 1:N, imag(ew), '.r')
	title 'mode frequencies', legend real imag Location Best
ew = real(ew);
figure, semilogy(1:N-1, ew(2:end)-ew(1:end-1), '+k')
	title 'mode frequency steps'
figure, imagesc(abs(X)), axis image, colormap gray, colorbar
title 'raw self products'
figure, imagesc(abs(r.dV*(U.'*V - V.'*U))), axis image, colormap gray, colorbar
title 'raw cross products'
figure, imagesc(abs(r.dV*(UU'*UU - VV'*VV))), axis image, colormap gray, colorbar
title 'corrected self products (should be Id)'
figure, imagesc(abs(r.dV*(UU.'*VV - VV.'*UU))), axis image, colormap gray, colorbar
title 'corrected cross products (should be zero)'

% checks
