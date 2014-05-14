# Standard imports for an environment like Matlab

# This doesn't import any interpolation routines: Python keeps reinventing that wheel, going from square to pentagonal to triangular.  See scipy.ndimage.interpolation and scipy.interpolate

from numpy import \
	abs, amax, amin, all, allclose, arange, argmax, array, asarray, \
	ceil, concatenate, diag, diagflat, diagonal, dot, empty, eye, \
	floor, indices, isnan, linspace, load, log, log10, logspace, \
	logical_not, logical_or, \
	mean, meshgrid, nan, ndarray, newaxis, nonzero, \
	pi, prod, ptp, \
	sign, sqrt, zeros, \
	exp, sin, cos, tan, tensordot
from numpy.random import rand, randn, seed
from numpy.fft import fft, ifft, fftn, ifftn, fftfreq
from numpy.linalg import eig, norm, qr

# work around Pyenv breaking the OSX backend for Matplotlib
import matplotlib; matplotlib.use("Agg")
from matplotlib.pyplot import colorbar, contour, contourf, figure, imshow, legend, loglog, plot, savefig, semilogy, subplot, text, title, xlabel, ylabel, xlim, ylim, xticks, yticks

seed(666)		# repeatability