function in = cogs(in)
%COGS     XPDE structure to sample a coherent order parameter
%    See BOGS for explanation.

out.name = 'coherent initial state';

% equilibrium order parameter by imaginary time

in.randoms = 2;
in.initial = @(w,r) reshape(repmat(r.a.op,1,r.ensembles(1)), r.d.a) + [1 1i]*w/2;

end
