function out = twop(in)
%TWOP     XPDE input structure for TW sampled coherent field states
% out.transfer adds Wigner noise to the input order parameter
% this assumes that the order parameter is normalised to the particle number

out.name = 'coherent initial state';
for s = {'dimension', 'fields', 'ranges', 'c', 'a', 'points'}
	if isfield(in, s{1}), out.(s{1}) = in.(s{1}); end
end

out.points(1) = 2;  out.ranges(1) = 1;
out.noises = 2;
out.transfer = @tfr;

end

function a = tfr(w,r,a0,r0)
	
	a = a0 + [1 1i]*w/2;
end