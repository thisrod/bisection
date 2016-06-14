% Find a bogoliubov ground state for the Vienna trap

op = static(vienna('trap'), 0);
op.ranges(1) = 0.5;
op.points(1) = 49;
op.steps = 30;
op = eqop(op);

uvs = bogs(op);

