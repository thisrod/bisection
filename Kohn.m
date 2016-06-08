function Kohn()

op = static(vienna('trap'), 0);
op.ranges(1) = 0.5;
op.points(1) = 49;
op.steps = 30;
op = eqop(op);

kohn = static(vienna('trap'), 0);
kohn.name = 'dynamics of z-axis Kohn mode';
kohn.ranges(1) = 50;
kohn.transfer = @(w,r,a0,r0) exp(1i*0.01*r.x).*a0;
kohn = xinstrument(kohn, 'n', 10, 'x', 'kx');
kohn.steps = 150;

xspde({cogs(op), kohn})

end
