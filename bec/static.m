function in = static(in, t)
%STATIC	Transform a dynamic atomic trap structure to be stuck at a fixed time
%
%    static(in, t)

in.name = sprintf('%s at time %f', in.name, t);
K = in.a.K;  in.a.K = @(r) K(setfield(r, 't', t));

end