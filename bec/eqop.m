function out = eqop(in)
%EQOP     annotate input structure with equilibrium order parameter
%
%    in = EQOP(in) annotates in.a.  This can be copied to other input structures.
%
%    See also: TRAP

out = in;
in = order(in);
in.raw = true;
in = xinstrument(in, 'n', 'N', 'Kx', 'rshl');

[~, in, rslt, raw] = xsim(in);  rslt = rslt{1};
fprintf('raw class: %s, size: %s\n', class(raw), mat2str(size(raw)))
fprintf('raw{1,2} class: %s, size: %s\n', class(raw{1,2}), mat2str(size(raw{1,2})))
a = squeeze(raw{1,2}(1,1,end,:));  a = a(:);  out.a.op = a;
kix = find(strcmp(in.olabels, 'K^2'));  out.a.Kludge = squeeze(rslt{kix}(1,1,end,:));
% Beware: xave fakes its answer to work around the grid length design flaw,
% so xave(uniform_value) is not equal to the uniform value!
out.a.healing = 1/sqrt(out.a.g*a(1)^2);

out.dV = in.dV;		% bogs needs this

end