function in = static(in, t)
%STATIC	Transform a dynamic atomic trap structure to be stuck at a fixed time

in.name = sprintf('%s at time %f', in.name, t);
J = in.c.K;
in.c.K = @(r) J(setfield(r, 't', t));

end