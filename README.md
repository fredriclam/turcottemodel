# turcottemodel
Finite volume solution for isothermal Euler of a mixed magma-vapour mixture model.

Based on [Turcotte et al. 1990](https://doi.org/10.1111/j.1365-246X.1990.tb01763.x), this learning code models mixture mass and momentum conservation laws, and conservation of the melt phase. The primitive variables are mixture density, momentum density using the mixture density, and volume fraction of melt (with dissolved gas).

Three basic finite volume methods are implemented in the code: (1) a piecewise constant reconstruction with flux computed in quasilinear form, upwinded using a matrix split of the flux jacobian; (2) a piecewise linear reconstruction with Lax-Friedrichs flux; and (3) a piecewise constant reconstruction with Local Lax-Friedrichs (Rusanov) flux.
