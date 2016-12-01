% Hartree orbitals in the trap plane, for fitted quartic potentials

% beware: qrtfit sets t

% Spectral matrix from Program 8 of 2000-Trefethen-Spectral

n1 = 7000/200;	% 1D density

% y axis samples N points in (-L,L], x axis 3N in (-3L,3L]
% keep the axes consistent with qrtfit.m
qrtfit,  L = 1;  N = 16;  h = 2*L/N;  
x = h*(1:3*N)-3*L;  y = h*(1:N)-L;  [X,Y] = ndgrid(x,y);

dxx = D2(3*N, 3*L);  Dxx = kron(eye(N), dxx);
dyy = D2(N, L);  Dyy = kron(dyy, eye(3*N));
LAP = Dxx + Dyy;

hit = 6;					% Hartree iterations
K = zeros(18, 3*N^2);
q = zeros(18, 3*N^2, 2);		% ground state and first excited state
w = zeros(18, hit, 3*N^2, 2);

for n = 1:2, for t = 0:17
	K(t+1,:) = [X(:).^4 ones(size(X(:))) X(:).^2 Y(:).^2]*cfs(1:4,t+1);
	for j = 2:hit
		HAM = 0.5*(-LAP+diag(K(t+1,:)));
		% average pairs of densities for stability
		% n2 = squeeze(abs(w(t+1, j-1, :, n)).^2 + abs(w(t+1, j-2, :, n)).^2) / 2;
		n2 = squeeze(abs(w(t+1, j-1, :, n)).^2);
		disp([n1 h^2*sum(n2(:)) max(n2(:)) n1*max(n2(:))])
		HAM = HAM + diag(0.5*0.1330*n1*n2);
		[Q, ks] = eig(HAM);
		[~,i] = sort(diag(ks));
		w(t+1, j, :, n) = Q(:,i(n))/h;
	end
end, end

% gnuplot can't do contours, so save for Matlab

save -mat horbsplotdat.mat t x y K q w
