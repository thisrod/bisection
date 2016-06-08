% tests of spectral differentiation

%% test that sin'' = -sin in 1D

n = 16;

in.dimension = 2;
in.points = [1 n];
in.ranges = [1 2*pi*(1-1/n)];
in.initial = @(~,r) sin(r.x);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
assert(norm(s+ssd(out)*s) < 1e-10);

%% test that sin(x) has laplacian -sin(x) in 2D

n = 16;

in.dimension = 3;
in.points = [1 n n];
in.ranges = [1 2*pi*(1-1/n) 1];
in.initial = @(~,r) sin(r.x);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
assert(norm(s+ssd(out,2)*s) < 1e-10);
