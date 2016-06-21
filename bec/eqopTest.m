clear notrap hmtrap box

notrap.name = 'no trap';
notrap.a.gamma = 1;
notrap.a.N = 3;
notrap.points = [49 70];
notrap.ranges = [10 69/70];		% work around the xint circular grid idiocy
notrap.steps = 30;
notrap = trap(notrap);

% watch out: the large potential, combined with XSPDE's invisible extrapolation, can give bizarre results when the fine grid has converged but the coarse one is unstable.

box.name = 'unit length box';
box.a.g = 0;
box.a.K = @(r) 1e3*(abs(r.x)>1/2);
box.a.N = 1;
box.points = [49 70];
box.ranges = [10 2];
box.steps = 60;
box = trap(box);

hmtrap.name = 'harmonic with l_0 = 1';
hmtrap.a.g = 0;
hmtrap.a.K = @(r) r.x.^2;
hmtrap.a.N = 1;
hmtrap.points = [49 70];
hmtrap.ranges = [10 20];
hmtrap.steps = 30;
hmtrap = trap(hmtrap);

% use debug option to eqop to get graphs

notrap = eqop(notrap);

hmtrap = eqop(hmtrap);

box = eqop(box);

%% repulsion energy for free atoms

assert(abs(notrap.a.Reqm - 18) < 1e-5);

%% potential energy is small in a box

assert(abs(box.a.Keqm - 0) < 10);

%% kinetic energy for the box ground state

assert(abs(box.a.Teqm - pi^2) < 1)

%% potential energy for harmonic

assert(abs(hmtrap.a.Keqm - 0.5) < 1e-4);

%% kinetic energy for harmonic

assert(abs(hmtrap.a.Teqm - 0.5) < 1e-4);
