
QuadSEM is specifically designed to compute the elastic full Hessian kernels on the fly based on seismic wave simulation using the spectral-element method (SEM). It derives its name from its capability to solve four elastic wave equations simultaneously based on the SEM. The development of QuadSEM is built upon the SPECFEM2D package (https://geodynamics.org/cig/software/specfem2d/, last accessed: 08.2019).

The SPECFEM2D is typically used with a single input model (m1), where each GLL point of the mesher stores one model value (e.g., Vp1, Vs1, and Rho1) for that point. In contrast, QuadSEM works with two models (m1 and m2), where each GLL point stores two model values instead of one. This allows the code to compute two sets of forward fields, s(m1) and s(m2), simultaneously during the forward simulation, and four wavefields simultaneously during the adjoint simulation: two (reversed or reconstructed) forward fields, s(m1) and s(m2), and two adjoint fields, s*(m1) and s*(m2). These wavefields are computed on the fly for each time step, enabling the computation of full Hessian kernels on the fly and within the same memory space. This approach avoids the need for large wavefield storage and extensive I/O operations for full Hessian kernel calculations.


Installation and Usage: 
 
1. Install the codes (tested only with 'gfortran gcc' or 'ifort and icc'):
./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi --with-scotch-dir=/home/QuadSEM_elastic/external_libs/scotch_5.1.12b  
or
./configure FC=ifort CC=icc MPIFC=mpif90 --with-mpi --with-scotch-dir=/home/QuadSEM_elastic/external_libs/scotch_5.1.12b
make clean
make all


 2. Using the Codes (It is better to be familiar with SPECFEM2D in advance)
     To simplify the usage, we have kept it almost the same as SPECFEM2D, where you only need to set three additional parameters in the Par_file for both forward and adjoint simulations (see the Par_file). The usage of the adjoint source is slightly different, as you need to add a suffix '_m2' for the adjoint source computed for m2. For example, AA.S****.BXX.adj for the first model and AA.S****.BXX.adj_m2 for the second model. The usage of other parameters remains the same as in SPECFEM2D. In the first version, we have implemented only the following: time_stepping_scheme = 1 and STACEY_ABSORBING_CONDITIONS = .true.
     To simplify the usage, we have adapted the same forward and adjoint simulation commands (e.g., mpirun -np 4 ./bin/xmeshfem2D; mpirun -np 4 ./bin/xspecfem2D) exactly as done in SPECFEM2D. However, please note that the forward and adjoint simulation in QuadSEM works for two models (m1 and m2) simultaneously when the associated QuadSEM parameters are set to 'true' in the Par_file. The Frechet and full Hessian kernels are computed during the adjoint simulation.


 3. Compute adjoint sources:
    Copy the seismogram of one model into the './OUTPUT_FILES_SEM' directory.
    Run './bin/xadj_seismogram t_min t_max AA.S0001 1' to obtain the x-component of the adjoint source, which will be stored in the './SEM/' directory.
    Run './bin/xadj_seismogram t_min t_max AA.S0001 3' to obtain the z-component of the adjoint source, which will also be stored in the './SEM/' directory. Please note that the x-component obtained in the previous step will be replaced by this computation.
The same steps apply to the seismograms of m2. It is important to add the suffix '_m2' for the adjoint source computed for m2, for example, 'AA.S****.BXX.adj_m2'


4. Try the above steps for two simple models as shown in the reference below, and you will understand the advantages of using this method.
    If you have any issues related to the implementation and usage, please contact the authors. More updates and examples will be provided in our next version. Occasionally, you may encounter a few invalid encodings in the *.jpg output, but this is only for display purposes and does not affect the wavefields and kernels. We have verified that the wavefields, waveforms and kernels are computed exactly the same as in SPECFEM2D when combining four solvers for two input models. The two input models are set to be identical to the external model used by SPECFEM2D (refer to Par_file).


------
Yujiang Xie
Ocean and Earth Science
University of Southampton, UK

Reference:
Yujiang Xie, Catherine A. Rychert, Nicholas Harmon, Qinya Liu and Dirk Gajewski. On‐the‐Fly Full Hessian Kernel Calculations Based upon Seismic‐Wave Simulations, Seismological Research Letters 2021,92,3832-3844, doi: https://doi.org/10.1785/0220200410







