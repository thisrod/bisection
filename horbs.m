# Hartree orbitals in the trap plane

# Spectral matrix from Program 8 of 2000-Trefethen-Spectral

prinax

[x, y, z, K] = loadk(0, U, r0);
[~,j0] = min(K, [], 2);  j0 = repmat(j0, [1,3,1]);
d = size(j0);
[i,j,k] = ndgrid(1:d(1), 1:d(2), 1:d(3));	% is there a special function for this?
K3 = K(sub2ind(size(K), i(:), j(:)+j0(:), k(:)));  K3 = reshape(K3,d);
K2 = K3(:,2,:) + (K3(:,1,:)-K3(:,3,:)).^2/8./(2*K3(:,2,:)-K3(:,1,:)-K3(:,3,:));
figure; contour(squeeze(K2), (0:0.1:1).^2, '-k')

if 1==0
# construct a spectral transverse laplacian operator
ns = size(K);  ns = ns([1 3]);  hs = 2*pi./ns;
Ls = [range(x) range(z)].*ns/2./(ns-1);
for i = 1:2
	column = [-pi^2/(3*hs(i)^2)-1/6 ...
		-0.5*(-1).^(1:ns(i)-1)./sin(hs(i)*(1:ns(i)-1)/2).^2];
	D2 = (pi/Ls(i))^2*toeplitz(column);
	switch i
		case 1
			dxx = D2;
			Dxx = kron(eye(ns(2)), D2);
		case 2
			dzz = D2;
			Dzz = kron(D2, eye(ns(1)));
	endswitch
endfor
LAP = Dxx + Dzz;
q = zeros(size(K(:)));
HAM = 0.5*(-LAP+diag(K(:)+0.283*abs(q).^2));
[Q, ks] = eig(HAM);
[k,i] = min(ks);

endif		% 1==0