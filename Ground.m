function Ground()                        %%name of main function

vienna obervables
order = vienna('initial');
order = xinstrument(order, 'T', 'K', 'V');

ground = bdv(order, 0);
ground = xinstrument(ground, 'Re', 2);

xspde({order, ground})

end
