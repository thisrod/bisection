% working script for Bogoliubov mode orthogonalisation

N = r.nspace -1;
keys = round(ew, 10);		% identify ews that are within 1e-10
L = r.ranges(2)*(r.points(2)+1)/r.points(2);	% circular length
r.a.healing = sqrt(L/r.a.Reqm);
kk = 2*pi*floor(1:0.5:r.points(2)/2)/L;	% bdg mode wave numbers

X = r.dV*(U'*U - V'*V);
Un = r.dV*sum(abs(U).^2);  Vn = r.dV*sum(abs(V).^2);

UU = nan(size(U));  VV = nan(size(V));

for key = unique(keys)'
	ixs = find(keys == key);
	if numel(ixs) == 1
		UU(:,ixs) = U(:,ixs);  VV(:,ixs) = V(:,ixs);
		continue
	end
	assert(numel(ixs) == 2);
	j = ixs(1);  k = ixs(2);
	numerator = -X(j,k) + sqrt(X(j,k)^2-(Un(j)-Vn(j))*(Un(k)-Vn(k)));
	l = numerator/(Un(k)-Vn(k));
	m = numerator/(Un(j)-Vn(j));

	UU(:,j) = U(:,j) + l*U(:,k);  VV(:,j) = V(:,j) + l*V(:,k);  
	UU(:,k) = U(:,k) + m*U(:,j);  VV(:,k) = V(:,k) + l*V(:,j);  
end

UUn = r.dV*sum(abs(UU).^2);  VVn = r.dV*sum(abs(VV).^2);

c = UUn - VVn;
c = 1./sqrt(c);
UU = UU.*repmat(c,r.nspace,1);  VV = VV.*repmat(c,r.nspace,1);  

% new_modes = new_modes.*repmat(c,r.nspace,2,1);
% 
% UU = squeeze(new_modes(:,1,:));  VV = squeeze(new_modes(:,2,:));

% Check that BM satisfies the equation it should
Prdl = (H0^2 + 2*r.a.g*diag(abs(a).^2))*H0*BP - BP*diag(ew.^2);
Mrdl = (H0^2 + 2*r.a.g*H0*diag(abs(a).^2))*BM - BM*diag(ew.^2);
figure, plot(1:N, sqrt(sum(abs(Prdl).^2)), '+k', ...
	1:N, sqrt(sum(abs(Mrdl).^2)), 'ok')
title 'residuals', legend '\psi^+' '\psi^-'

figure, plot(1:N, real(ew), 'ok', 1:N, imag(ew), '.r', 1:N, sqrt(2)*kk/r.a.healing, '+k')
	title 'mode frequencies', legend real imag 'dspn relation' Location Best

% confirm that orthogonalised (u,v) satisfy the BdG problem?

Urhs = H0*UU+r.a.g*diag(abs(a).^2)*(UU-VV);
Vrhs = H0*VV+r.a.g*diag(abs(a).^2)*(VV-UU);

figure, semilogy(1:N, sqrt(sum(abs(Urhs-UU*diag(ew)))) ./ sqrt(UUn+VVn), '^k', ...
	1:N, sqrt(sum(abs(Vrhs+VV*diag(ew)))) ./ sqrt(UUn+VVn), 'vk')
title 'fractional residuals of BdG problem', legend u v

figure, imagesc(abs(r.dV*(U'*U - V'*V))), axis image, colormap gray, colorbar
title 'raw self products'

figure, imagesc(abs(r.dV*(U.'*V - V.'*U))), axis image, colormap gray, colorbar
title 'raw cross products'

figure, imagesc(abs(r.dV*(UU'*UU - VV'*VV))), axis image, colormap gray, colorbar
title 'corrected self products (should be Id)'

figure, imagesc(abs(r.dV*(UU.'*VV - VV.'*UU))), axis image, colormap gray, colorbar
title 'corrected cross products (should be zero)'

figure, plot(1:N, nms(1,:), 'ok', 1:N, nms(2,:), 'xk', [1:N; 1:N], nms, '-k' );
	title 'norms of U and V modes'
figure, plot(1:N, real(ew), '.k', 1:N, imag(ew), '.r')
	title 'mode frequencies', legend real imag Location Best
ew = real(ew);
figure, semilogy(1:N-1, ew(2:end)-ew(1:end-1), '+k')
	title 'mode frequency steps'

% checks
