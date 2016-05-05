function in = order(in, N, t)
%ORDER	Transform an atomic trap structure to find its ground state order parameter
%
%    out = ORDER(in, N, t) computes the order parameter for N atoms in the trap potential at time t
%
%    out = ORDER(in, N) assumes the trap potential is constant
%
% The user will add the duration of the imaginary time integration later as in.ranges(1).

if nargin == 3, in = static(in, t); end
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
