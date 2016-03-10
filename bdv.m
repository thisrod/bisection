function out = bdv(in)
%BDV     XPDE input structure for Bogoliubov ground states
% For now, set the field to the lagrangian of the input field.

out.name = 'Bogoliubov initial state';
for s = {'dimension', 'fields', 'ranges', 'c', 'points'}
	out.(s{1}) = in.(s{1});
end

out.points(1) = 2;  out.ranges(1) = 1;
out.transfer = @tfr;

end

function a = tfr(w,r,a0,r0)
	size(a0)
	p = r.points,  g = r.ranges
	Dxx = kron(eye(p(3)), ssd(p(2), g(2)));
	Dyy = kron(ssd(p(3), g(3)), eye(p(2)));
	a = (Dxx+Dyy)*a0.';
	a = a.';
end