# PDEFourierAnalysis.jl

This is a package for calculating and plotting the stability regions of
numerical methods via Fourier analysis. The aim is to provide many of the
standard quantities of interest in the numerical analysis for spatial,
temporal and coupled methods.

# Organisation
Methods are broken down into three groups, these are:
 - Temporal, which are ODE methods
 - Spatial, which are methods that approximate spatial derivatives of a given PDE
 - Scheme, which are the combination of a spatial and temporal method
