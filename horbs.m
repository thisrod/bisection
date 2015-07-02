% Hartree orbitals in the trap plane, for fitted quartic potentials

% TODO check 2D repulsion aa
% beware: qrtfit sets t

% Spectral matrix from Program 8 of 2000-Trefethen-Spectral

% y axis samples N points in (-L,L], x axis 3N in (-3L,3L]
aa = 0.1330*7000/200;
qrtfit,  t = 17;  L = 1;  N = 16;  h = 2*L/N;  
x = h*(1:3*N)-3*L;  y = h*(1:N)-L;  [X,Y] = ndgrid(x,y);

K = [X(:).^4 ones(size(X(:))) X(:).^2 Y(:).^2]*cfs(1:4,t+1);

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

q = zeros(3*N^2, 5);
for l = 2:5
	HAM = 0.5*(-LAP+diag(K+0.1*l*aa*abs(q(:,l-1)).^2));
	[Q, ks] = eig(HAM);
	[k,i] = min(diag(ks));
	q(:,l) = sqrt(7000)*Q(:,i);
endfor

% gnuplot can't do contours, so save for Matlab

ql = reshape(q(:,2), size(X));
save -mat horbsplotdat.mat t x y K ql
