# Multi-SEM
DOI:10.5281/zenodo.4075546

codes are available from the author, to be uploaded when ready for acceptance!

In this Multi-SEM, we implement 4 methods to compute the Hessian kerenls:

1. The on-the-fly implementation (for elastic case only):
2. The parsimonious storage method (for elastic and anelastic case):
3. The full wavefield storage method (for elastic and anelastic case):
4. A finite-difference approximation (for elastic and anelastic case): 

With importing two models (GLL or ASCII format), the Multi-SEM can compute these kernels simultaneously with only requiring about a 2~2.5 computational cost required in comparison to the computation of Frechet kernel alone:
1) the Frechet kernels for m1 and m2
2) the finite-difference approximation based on g(m1) and g(m2)
3) the approximate Hessian kernel (to be used with Gaussian-Newton FWI)
4) the full Hessian kernel (to be used with Newton FWI)

The Multi-SEM currently supports computing the Hessian kernel using two models, and the extension for models > 3 is straightforward, where one can compute the forward and adjoint fields for multiple stets of models simultaneously. 

For the wave simulation, we use the spectral-element and adjoint methods (Tromp et al., 2005). The same ideas can be used for other wave-equation solvers as well. 

----Yujiang Xie, 2021. 

