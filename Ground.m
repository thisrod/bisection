function Ground()                        %%name of main function

vienna obervables
order = vienna('initial');

ground = bdv(order);
ground = xinstrument(ground, 'n', 2);

xspde({order, ground})

end
