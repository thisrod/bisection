function D2=ssd(in, axis)
%SSD     Spectral differentiation matrix for XPDE input structure
%
%    D2 = SSD(in, axis) sets D2 to a matrix that, when applied to a field
%    of shape in.d.a, gives the second spatial derivative of that field along the numbered axis.
%    This is useful for explicit diagonalisation and so on.
%
%    D = SSD(in, 'lap') returns the laplacian.
%
%    The axis argument can be omitted when the field has only one spatial dimension.
%
%    See Trefethen, "Spectral Methods in Matlab".
	
	if nargin == 1 && in.dimension == 2
		axis = 2;
	elseif nargin == 1
		error 'No axis specified along which to take the derivative'
	end

	% Laplacian by recursion
	if strcmp(axis, 'lap')
		D2 = 0;
		for i = 2:in.dimension, D2 = D2 + ssd(in,i); end
		return
	end	

	% calculate second derivative operator for chosen axis
	N = in.points(axis);  L = in.ranges(axis)/2;
	assert(~mod(N,2));	% odd grids NYI
	L = L*N/(N-1);  h = 2*pi./N;
	column = [-pi^2/(3*h^2)-1/6 ...
		-0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
	D1 = (pi/L)^2*toeplitz(column);

	% extend to an operator over all axes
	
	D2 = 1;
	for i = 2:in.dimension
		if i == axis
			D2 = kron(D1, D2);
		else
			D2 = kron(eye(in.points(i)), D2);
		end
	end

end
