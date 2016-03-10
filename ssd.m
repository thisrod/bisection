function D2=ssd(N,L)
	% spectral second derivative, N points on (-L,L)
	% note changed to match XSPDE grids
	L = L*N/(N-1);  h = 2*pi./N;
	column = [-pi^2/(3*h^2)-1/6 ...
		-0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
	D2 = (pi/L)^2*toeplitz(column);
end
