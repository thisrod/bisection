% working script for Bogoliubov mode orthogonalisation

N = r.nspace -1;
L = r.ranges(2)*(r.points(2)+1)/r.points(2);	% circular length
r.a.healing = sqrt(L/r.a.Reqm);
kk = 2*pi*floor(1:0.5:r.points(2)/2)/L;	% bdg mode wave numbers
U = W(1:end/2,:);  V = W(end/2+1:end, :);
F = r.dV*(U'*U - V'*V);  G = r.dV*(U'*V-V'*U);
FF = r.dV*(UU'*UU - VV'*VV);  GG = r.dV*(UU'*VV-VV'*UU);
oix = find(max(abs(GG))>1);
BPn = sqrt(r.dV*sum(abs(BP).^2));  BMn = sqrt(r.dV*sum(abs(BM).^2));
Un = r.dV*sum(abs(U).^2);  Vn = r.dV*sum(abs(V).^2);
UUn = r.dV*sum(abs(UU).^2);  VVn = r.dV*sum(abs(VV).^2);


% Check that BM satisfies the equation it should
Prdl = (H0^2 + 2*r.a.g*diag(abs(a).^2)*H0)*BP - BP*diag(ew.^2);
Mrdl = (H0^2 + 2*r.a.g*H0*diag(abs(a).^2))*BM - BM*diag(ew.^2);
figure, plot(1:N, sqrt(sum(abs(Prdl).^2)), '+k', ...
	1:N, sqrt(sum(abs(Mrdl).^2)), 'ok')
title 'residuals', legend '\psi^+' '\psi^-'

figure, plot(1:N, real(ew), 'ok', 1:N, imag(ew), '.r', 1:N, sqrt(2)*kk/r.a.healing, '+k')
	title 'mode frequencies', legend real imag 'dspn relation' Location Best

% pattern of imaginary mode functions

imixs = find(~all(neareal(V)));
figure, subplot(length(imixs)+2, 1,1), spy(~all(neareal(U)))
	title 'Complex mode functions'
subplot(length(imixs)+2, 1,2), spy(~all(neareal(V)))
for i = 1:length(imixs)
	j = imixs(i);
	subplot(length(imixs)+2, 1,i+2)
	plot(r.x, real(U(:,j)), '-k', r.x, imag(U(:,j)), '-r', ...
		r.x, real(V(:,j)), '--k', r.x, imag(V(:,j)), '--r')
end

% confirm that orthogonalised (u,v) satisfy the BdG problem?

Urhs = H0*UU+r.a.g*diag(abs(a).^2)*(UU-VV);
Vrhs = H0*VV+r.a.g*diag(abs(a).^2)*(VV-UU);

figure, semilogy(1:N, sqrt(sum(abs(Urhs-UU*diag(ew)))) ./ sqrt(UUn+VVn), '^k', ...
	1:N, sqrt(sum(abs(Vrhs+VV*diag(ew)))) ./ sqrt(UUn+VVn), 'vk')
title 'fractional residuals of BdG problem', legend u v

figure
subplot 221
imagesc(abs(F)), axis image, colormap gray, colorbar
title 'raw self products'

subplot 222
imagesc(abs(G)), axis image, colormap gray, colorbar
title 'raw cross products'

subplot 223
imagesc(abs(FF)), axis image, colormap gray, colorbar
title 'corrected self products (should be Id)'


subplot 224
imagesc(abs(GG)), axis image, colormap gray, colorbar
title 'corrected cross products (should be zero)'

figure, plot(1:N, sqrt(UUn), 'ok', 1:N, sqrt(VVn), 'xk');
	title 'norms of corrected U and V modes', legend U V

ew = real(ew);
figure, semilogy(1:N-1, ew(2:end)-ew(1:end-1), '+k')
	title 'mode frequency steps'

% checks

figure, semilogy(1:N, BPn, '^k', 1:N, BMn, 'vk')
	title 'norms of \psi^+ and \psi^-'

figure
subplot 221, plot(r.x, U(:,oix(1)), '-k', r.x, V(:,oix(1)), '-r')
	title 'original mode a'
subplot 222, plot(r.x, U(:,oix(2)), '-k', r.x, V(:,oix(2)), '-r')
	title 'original mode b'
subplot 223, plot(r.x, UU(:,oix(1)), '-k', r.x, VV(:,oix(1)), '-r')
	title 'corrected mode a'
subplot 224, plot(r.x, UU(:,oix(2)), '-k', r.x, VV(:,oix(2)), '-r')
	title 'corrected mode b'
