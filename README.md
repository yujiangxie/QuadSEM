# QuadSEM

Quad Spectral-Element Method (QuadSEM) is a recently developed solver designed to simultaneously solve four elastic wave equations for two models, where each wave equation is solved by the classical SEM using the SPECFEM2D package (https://geodynamics.org/cig/software/specfem2d/). The QuadSEM is initially designed to compute the full Hessian kernels or full Hessian vector products on the fly using the exact forward and adjoint fields on the fly, substantially reducing the disk space and I/O requirements for storing and transmitting the forward and adjoint fields for the full Hessian kernel calculations due to the on-the-fly feature. The QuadSEM can be used to simultaneously compute the Frechet and full Hessian kernels requiring only about a 2-fold computational cost in comparison to the computation of the Frechet kernels alone.

The QuadSEM is a development of the Specfem2D, accounting for two input models m1 and m2, so that people can compute two sets of forward fields s(m1) and s(m2) simultaneously in the forward simulation, as well as can compute four sets of wave fields simultaneously in the adjoint simulation: that is two (reversed) forward fields s(m1) and s(m2) and two adjoint fields s*(m1) and s*(m2). We compute these fields simultaneously on the fly so that the full Hessian kernels can be computed on the fly as well since the construction of the full Hessian kernels needs these fields. 

For better usage, we use the framework of the Specfem2D that many people are familiar with, where the user just needs to set four parameters in the Par_file for the forward and adjoint simulation. The usage of the other parameters is the same as that of the Specfem2D.  

All codes are given in QuadSEM2D_elastic.zip
