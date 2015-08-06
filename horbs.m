% Hartree orbitals in the trap plane, for fitted quartic potentials

% beware: qrtfit sets t

% Spectral matrix from Program 8 of 2000-Trefethen-Spectral

q1 = sqrt(7000/200);	% 1D order parameter

% y axis samples N points in (-L,L], x axis 3N in (-3L,3L]
% keep the axes consistent with qrtfit.m
qrtfit,  L = 1;  N = 16;  h = 2*L/N;  
x = h*(1:3*N)-3*L;  y = h*(1:N)-L;  [X,Y] = ndgrid(x,y);

function D2=ssd(N,L)
	% spectral second derivative, N points on (-L,L]
	h = 2*pi./N;
	column = [-pi^2/(3*h^2)-1/6 ...
		-0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
	D2 = (pi/L)^2*toeplitz(column);
endfunction

dxx = ssd(3*N, 3*L);  Dxx = kron(eye(N), dxx);
dyy = ssd(N, L);  Dyy = kron(dyy, eye(3*N));
LAP = Dxx + Dyy;

hit = 5;					% Hartree iterations
K = zeros(17, 3*N^2);
q = zeros(17, 3*N^2, 2);		% ground state and first excited state
w = zeros(17, hit, 3*N^2);

for t = 0:17
	K(t+1,:) = [X(:).^4 ones(size(X(:))) X(:).^2 Y(:).^2]*cfs(1:4,t+1);
	for j = 2:hit
		HAM = 0.5*(-LAP+diag(K(t+1,:)));
		ww = squeeze(w(t+1, j-1, :));
		HAM = HAM + diag(0.5*0.1330*q1*abs(ww).^2);
		[Q, ks] = eig(HAM);
		[~,i] = sort(diag(ks));
		w(t+1, j, :) = sqrt(7000)*Q(:,i(1));
	end
end

% gnuplot can't do contours, so save for Matlab

save -mat horbsplotdat.mat t x y K q w
