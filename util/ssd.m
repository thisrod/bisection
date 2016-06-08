function D2=ssd(in)
%SSD     Spectral differentiation matrix for XPDE input structure
%
%    D2 = SSD(in) sets D2 to a matrix that, when applied to a field
%    of shape in.d.a, gives the second spatial derivative of that field.
%    This is useful for explicit diagonalisation and so on.

	N = in.points(2);  L = in.ranges(2)/2;
	L = L*N/(N-1);  h = 2*pi./N;
	column = [-pi^2/(3*h^2)-1/6 ...
		-0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
	D2 = (pi/L)^2*toeplitz(column);
end
