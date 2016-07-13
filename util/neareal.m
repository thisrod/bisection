function x = neareal(A, rtol)
%NEAREAL elements that are nearly real
%
%    Ask Rodney to fill this in.

if nargin < 2, rtol = 1e-5; end
x = abs(imag(A))/norm(A(:), inf)  > rtol;

end
