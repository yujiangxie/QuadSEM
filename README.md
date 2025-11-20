====== QuadSEM ======

The Quad Spectral-Element Method (QuadSEM) is primarily designed to compute the full Hessian kernels on the fly using the exact forward and adjoint fields. This approach significantly reduces the disk space and I/O requirements, often by more than three orders of magnitude or more depending on the time steps. By leveraging the on-the-fly feature, storing and transmitting the forward and adjoint fields for full Hessian kernel calculations is largely avoided. QuadSEM enables the simultaneous computation of both the Frechet and full Hessian kernels on the fly, with a computational cost approximately twice that of computing the Frechet kernels alone.

Hessian kernels have two primary applications: In inversion, they can enhance the convergence rate and mitigate trade-offs between inverted multi-parameters. After inversion, they can be utilized to analyze the resolution of the inverted model.

Compared to conventional wave equation solvers like SPECFEM2D, QuadSEM is specifically designed to accommodate two input models, namely m1 and m2. This allows for the simultaneous computation of two sets of forward fields, s(m1) and s(m2), during the forward simulation, as well as four sets of wave fields during the adjoint simulation. These include two reversed or reconstructed forward fields, s(m1) and s(m2), and two adjoint fields, s*(m1) and s*(m2). These fields are computed on the fly, eliminating the need for extensive wavefield storage. In the elastic case, only the boundary fields and the last snapshot of the forward fields are required.

In contrast to SPECFEM2D,, which employs a single model (e.g., vp, vs, and rho values) for simulations, QuadSEM utilizes two models (e.g., vp1, vp2, vs1, vs2, rho1, rho2) at each GLL point. This approach significantly reduces the intensive I/O communications between these fields compared to running SPECFEM2D twice in parallel.

While QuadSEM can be implemented in different computer languages or developed as a standalone software, we have chosen to leverage the SPECFEM2D framework for user convenience. Users only need to set four parameters in the Par_file for the new forward and adjoint simulations. The usage of other input parameters remains consistent with SPECFEM2D.

For the elastic codes, please also refer to QuadSEM2D_elastic.zip in the master branch

====== QuadSEM-Q ======

Feel free to contact the authors for access to the codes. Anyone who emails the authors can obtain the codes. We are currently working on a project, so the codes will be uploaded once the manual is ready.

Compared to SPECFEM2D, which employs a single model (e.g., vp, vs, rho, Qk, and Qu values) for each GLL point in the simulations, QuadSEM-Q utilizes two models (e.g., vp1, vp2, vs1, vs2, rho1, rho2, Qk1, Qk2, Qu1, Qu2) at each GLL point. This approach significantly reduces the need for extensive wavefield storage and minimizes the costly I/O communications between these fields, in contrast to running SPECFEM2D twice in parallel.

References:

Yujiang Xie, Catherine A. Rychert and Nicholas Harmon. Elastic and anelastic adjoint tomography with and full Hessian kernels, Geophysical Journal International, 234, 1205-1235, 2023, https://doi.org/10.1093/gji/ggad114
Yujiang Xie, Catherine A. Rychert, Nicholas Harmon, Qinya Liu and Dirk Gajewski. On‐the‐Fly Full Hessian Kernel Calculations Based upon Seismic‐Wave Simulations, Seismological Research Letters 2021,92, 3832-3844, doi: https://doi.org/10.1785/0220200410
