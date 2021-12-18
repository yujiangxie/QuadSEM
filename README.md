# ====== QuadSEM ======

The Quad Spectral-Element Method (QuadSEM) is initially designed to compute the full Hessian kernels or full Hessian vector products on the fly using the exact forward and adjoint fields, which substantially reduces the disk space and I/O requirements for storing and transmitting the forward and adjoint fields for the full Hessian kernel calculations due to the on-the-fly feature. The QuadSEM can be used to simultaneously compute the Frechet and full Hessian kernels on the fly requiring only about a 2-fold computational cost in comparison to the computation of the Frechet kernels alone.

Comparison to conventional wave equation solvers, the QuadSEM is designed for two input models m1 and m2, so that people can compute two sets of forward fields s(m1) and s(m2) simultaneously in the forward simulation, as well as can compute four sets of wave fields simultaneously in the adjoint simulation: that is two (reversed) forward fields s(m1) and s(m2) and two adjoint fields s*(m1) and s*(m2). These fields are simultaneously computed on the fly so that the full Hessian kernels can be computed on the fly as well without the need for the entire wavefield storage. Only the boundary fields and the last snapshot of the forward fields are needed in the elastic case.

The QuadSEM can be written in a different computer language or can be designed as new software. However, for easy using the codes, we use the framework of the Specfem2D that many people are familiar with, where one just needs to set four parameters in the Par_file for the new forward and new adjoint simulation. The use of the other input parameters is the same as that of the Specfem2D. 

All codes are given in QuadSEM2D_elastic.zip.
~~~~~~
Tips: compared to the Specfem2D, which uses one model (e.g., vp, vs, and rho) for the simulations, while in the QuadSEM, it uses two models (e.g., vp1,vp2, vs1,vs2, rho1,rho2, that is each GLL point has two values, instead of one).
~~~~~~

Reference:
Yujiang Xie, Catherine A. Rychert, Nicholas Harmon, Qinya Liu, Dirk Gajewski; On‐the‐Fly Full Hessian Kernel Calculations Based upon Seismic‐Wave Simulations. Seismological Research Letters 2021; doi: https://doi.org/10.1785/0220200410

# ====== QuadSEM-Q ======
The visco-elastic version will be uploaded when the paper is ready for acceptance. 
~~~~~~
Tips: compared to the Specfem2D, which uses one model (e.g., vp, vs, rho, Qk, Qu) for the simulations, while in the QuadSEM, it uses two models (e.g., vp1,vp2, vs1,vs2, rho1,rho2, Qk1,Qk2, Qu1,Qu2, that is each GLL point has two values, instead of one).
~~~~~~
