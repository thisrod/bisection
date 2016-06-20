% not yet set up as test cases

clear notrap hmtrap

notrap.name = 'no trap';
notrap.a.gamma = 1;
notrap.a.N = 1;
notrap.points = [49 70];
notrap.ranges = [10 nan];
notrap.steps = 30;
notrap = trap(notrap);

% find out why n goes negative

box.name = 'box';
box.a.g = 0;
box.a.K = @(r) 1e3*(abs(r.x)>1/2);
box.a.N = 1;
box.points = [49 70];
box.ranges = [10 2];
box.steps = 60;
box = trap(box);

hmtrap.name = 'harmonic';
hmtrap.a.g = 0;
hmtrap.a.K = @(r) r.x.^2;
hmtrap.a.N = 1;
hmtrap.points = [49 70];
hmtrap.ranges = [10 5];
hmtrap.steps = 30;
hmtrap = trap(hmtrap);

notrap = eqop(notrap, 'debug');

hmtrap = eqop(hmtrap, 'debug');

box = eqop(box, 'debug');

figure, plot(1:70, real(box.a.op), '-k', 1:70, imag(box.a.op), '-r')
