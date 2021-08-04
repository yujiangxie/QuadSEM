# QuadSEM

Quad Spectral-Element Method (QuadSEM) is a recently developed solver designed to simultaneously solve four elastic wave equations for two models, where each wave equation is solved by the classical SEM using the SPECFEM2D package (https://geodynamics.org/cig/software/specfem2d/). The QuadSEM is initially designed to compute the full Hessian kernels or full Hessian vector products on the fly using the exact forward and adjoint fields, substantially reducing the disk space and I/O requirements for storing and transmitting the forward and adjoint fields due to the on-the-fly feature. The QuadSEM can be used to simultaneously compute the Frechet and full Hessian kernels requiring only about a 2-fold computational cost in comparison to the computation of the Frechet kernels alone.


We are preparing the .doc files for the QuadSEM. Models and codes will be uploaded once the paper is available online....

