Programs in this directory

|prinax.m| determines the principle axes of the sampled potentials supplied from Vienna.  The trap axes don't align with the sampling axes.  It defines a function |loadk|, that loads the Vienna potential from the data, and |tcent|, that finds the approximate centre of the potential.

|qrtfit.m| generates the coefficients of quartic polynomials that fit the trap potentials.

Currently, the plotting facilities in Matlab work better than those in Octave.  In particular, Matlab can overlay contours on a zplot.  So the pattern is for qrtfit (say) to save data to be plotted in qrtfitdat.mat, which is loaded at the start of qrtfitfigs.m.