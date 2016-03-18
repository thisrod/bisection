% one dimensional BEC pair correlations

clear free
free.name = 'free bosons in one dimension';
free.a.gamma = 1e3;
free.a.N = 1e5;
free = trap(free);

ground = order(free, 1e5);
ground.ranges(1) = 0.5;
ground.points = [49 70];
ground.steps = 30;
ground = xinstrument(ground, 'n', 'N', 1e5);
xspde(ground)