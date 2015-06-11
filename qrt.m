# quartic potential, as suggested by Schiedmayer

L = 30;  N = 30;  h = 2*pi/N;  q = h*(1:N);  q = L*(q-pi)/pi;
K = 0.1*(q/10).^4 - 0.45*(q/10).^2;	# potential after bisection
D2K = 0.012*(q/10).^2 - 0.009;

column = [-pi^2/(3*h^2)-1/6 ...
	-0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
D2 = (pi/L)^2*toeplitz(column);

figure
subplot(211)
plot(q,K)
subplot(212)
plot(q,D2*K', q, D2K)
