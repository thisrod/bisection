% one dimensional BEC pair correlations

clear system gsop

system.name = 'free gas in one dimension';
system.a.gamma = 1e3;
system.a.N = 1e5;
system = trap(system);

gsop = order(system, 1e5);
gsop.ranges(1) = 0.5;
gsop.points = [49 70];
gsop.steps = 30;
gsop = xinstrument(gsop, 'n', 'N', 'T', 'K', 'V');

ground = system;
ground.ranges(1) = 0.5;
ground.points = [49 70];
ground = bdg(ground, 0);


xspde({gsop, ground})