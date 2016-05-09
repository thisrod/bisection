function in = order(in)
%ORDER	Transform an atomic trap structure to find its ground state order parameter
%
%    out = ORDER(in) computes the equilibrium order parameter the trap potential of in.  This should be a static potential, and in.a.N must be set to the number of atoms.
%
% The user must add the duration of the imaginary time integration as in.ranges(1).  This can be done before or after the call to order.

N = in.a.N;
in.name =	sprintf('Ground state order parameter with %d atoms for %s', N, in.name);
in.step = @(a,xi,r) nrmstp(a,xi,r,N);
in.initial = @(w,r) ones(size(r.x));
f = in.da;  g = in.linear;
in.da = @(a,w,r) -1i*f(a,w,r);
in.linear = @(D,r) -1i*g(D,r);

end

function b = nrmstp(a,xi,r,N)
	% normalise the order parameter for imaginary time integration
	b = xMP(a,xi,r);
	s = xint(abs(b).^2, r.dx, r);
	b = sqrt(N./s).*b;
end
