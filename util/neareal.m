function x = neareal(A, rtol)
%NEAREAL elements that are nearly real
%
%    NEAREAL(X) is a predicate for elements of X that are nearly real.
%    NEAREAL(X, RTOL) sets a limit on the imaginary part as a fraction of the maximum absolute value.


if nargin < 2, rtol = 1e-5; end
x = abs(imag(A))/norm(A(:), inf)  < rtol;

end
