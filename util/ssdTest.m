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

%% test that sin(x) has second derivative -sin(x) in 2D

n = 16;

in.dimension = 3;
in.points = [1 n n];
in.ranges = [1 2*pi*(1-1/n) 1];
in.initial = @(~,r) sin(r.x);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
assert(norm(s+ssd(out,2)*s) < 1e-10);

%% test that sin(y) has second derivative -sin(y) in 2D

n = 16;

in.dimension = 3;
in.points = [1 n n];
in.ranges = [1 1 2*pi*(1-1/n)];
in.initial = @(~,r) sin(r.y);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
assert(norm(s+ssd(out,3)*s) < 1e-10);

%% test that sin(3x+4y) has second derivative wrtx -9sin in 2D

n = 16;

in.dimension = 3;
in.points = [1 n n];
in.ranges = [1 2*pi*(1-1/n) 2*pi*(1-1/n)];
in.initial = @(~,r) sin(3*r.x+4*r.y);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
assert(norm(9*s+ssd(out,2)*s) < 1e-10);

%% test that sin(3x+4y) has second derivative wrty -16sin in 2D

n = 16;

in.dimension = 3;
in.points = [1 n n];
in.ranges = [1 2*pi*(1-1/n) 2*pi*(1-1/n)];
in.initial = @(~,r) sin(3*r.x+4*r.y);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
assert(norm(16*s+ssd(out,3)*s) < 1e-10);

%% test that sin(3x+4y) has laplacian -25sin in 2D

n = 16;

in.dimension = 3;
in.points = [1 n n];
in.ranges = [1 2*pi*(1-1/n) 2*pi*(1-1/n)];
in.initial = @(~,r) sin(3*r.x+4*r.y);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
LAP = ssd(out,'lap');
assert(norm(25*s+LAP*s) < 1e-10);

%% repeat with different number of points on each axis

m = 16;  n = 20;

in.dimension = 3;
in.points = [1 m n];
in.ranges = [1 2*pi*(1-1/m) 2*pi*(1-1/n)];
in.initial = @(~,r) sin(3*r.x+4*r.y);
in.raw = true;

[~,out,~,raw] = xsim(in);
s = squeeze(raw{1}(1,1,end,:));  s = s(:);
LAP = ssd(out,'lap');
assert(norm(25*s+LAP*s) < 1e-10);
