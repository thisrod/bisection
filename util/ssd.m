function D2=ssd(in, axis)
%SSD     Spectral differentiation matrix for XPDE input structure
%
%    D2 = SSD(in, axis) sets D2 to a matrix that, when applied to a field
%    of shape in.d.a, gives the second spatial derivative of that field along the numbered axis.
%    This is useful for explicit diagonalisation and so on.
%
%    The axis argument can be omitted when the field has only one spatial dimension.
%
%    See Trefethen, "Spectral Methods in Matlab".
	
	if nargin == 1 && in.dimension == 2
		axis = 2;
	elseif nargin == 1
		error 'No axis specified along which to take the derivative'
	end

	N = in.points(axis);  L = in.ranges(axis)/2;
	assert(~mod(N,2));	% odd grids NYI
	L = L*N/(N-1);  h = 2*pi./N;
	column = [-pi^2/(3*h^2)-1/6 ...
		-0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
	D2 = (pi/L)^2*toeplitz(column);
end
