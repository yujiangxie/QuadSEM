!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!=====================================================================

  subroutine compute_kernels() ! lucas, in the time loop: do it = 1,NSTEP 

! computes adjoint sensitivity kernel contributions
!
! see e.g. Tromp et al. (2005) for elastic calculation
! and Morency et al. (2009) for poroelastic calculation

  use constants, only: APPROXIMATE_HESS_KL

  use specfem_par, only: any_acoustic,any_elastic,any_poroelastic,CTD_SEM, Full_Hessian_by_Wavefield_Stored !lucas, CTD-SEM

  implicit none

  ! acoustic simulations
  if (any_acoustic) then
    call compute_kernels_ac()
  endif

  ! elastic simulations
  if (any_elastic) then
    
    if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas, CTD-SEM-------
     call compute_kernels_el_Ha_Hb_Hc_Habc() ! for elastic case
   !call compute_kernels_el_m2() ! for test only, not used
    else!---------------------------------
    call compute_kernels_el() !lucas, need to incoparate into below
    endif
  endif

  ! poro-elastic simulations
  if (any_poroelastic) then
    call compute_kernels_po()
  endif

  ! computes an approximative Hessian for preconditioning kernels
  if (APPROXIMATE_HESS_KL) then
    call compute_kernels_Hessian()
  endif

  end subroutine compute_kernels

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_el() !lucas, the old code and elastic
!lucas, for simulation_type=3, the displ_elastic is obtained when the source is at the receiver, and b_displ_elastic read from last frame of the forward
!lucas, the first it here is the last step of the forward.
! elastic kernel calculations
! see e.g. Tromp et al. (2005)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,HALF,TWO,FOUR_THIRDS

  use specfem_par, only: ispec_is_elastic,rho_k, & ! AXISYM,
                         rho_kl,mu_kl,kappa_kl,rhop_kl,beta_kl,alpha_kl,bulk_c_kl,bulk_beta_kl, &
                         nglob,nspec,ibool,b_displ_elastic,b_accel_elastic, & !lucas Ha
                         density,poroelastcoef,kmato,assign_external_model,rhoext,vsext,vpext, &
                         deltat,P_SV,displ_elastic, &
                         mu_k,kappa_k,ibool,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         GPU_MODE,it,NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS, &
                         islice_selected_rec, ispec_selected_rec, ibool,myrank ! lucas

                         

  use specfem_par_gpu, only: Mesh_pointer,deltatf

  use specfem_par, only: c11_k,c13_k,c15_k,c33_k,c35_k,c55_k,ispec_is_anisotropic, &
                         rho_kl,c11_kl,c13_kl,c15_kl,c33_kl,c35_kl,c55_kl

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob
  integer :: iglob_source_lucas
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal
!  real(kind=CUSTOM_REAL) :: rndx,rndy1,rndy2 !lucas


  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl
!  double precision :: rndmin, rndmax      ! lucas
!  rndmin = -0.0001
!  rndmax =  0.0001


  ! 1.lucas, to get kappa_k, mu_k, and rho_k at each iglob
  ! elastic kernels
  if (.not. GPU_MODE) then
    ! updates kernels on CPU
    do ispec = 1,nspec
      if (ispec_is_elastic(ispec)) then
        do j = 1,NGLLZ; do i = 1,NGLLX
          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL
          duz_dxi = 0._CUSTOM_REAL
          dux_dgamma = 0._CUSTOM_REAL
          duz_dgamma = 0._CUSTOM_REAL
          b_dux_dxi = 0._CUSTOM_REAL
          b_duz_dxi = 0._CUSTOM_REAL
          b_dux_dgamma = 0._CUSTOM_REAL
          b_duz_dgamma = 0._CUSTOM_REAL

          ! 1.1 lucas, to get dux_dxi, duz_dxi, dux_dgamma, duz_dgamma and also the back conterparts
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX 
            !u*(m1)
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
            !u(m1) 
            b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            b_duz_dxi = b_duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma = b_duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

!                !------- adding noise to u* and u for mimicing the compression method
!                call random_number(rndx)
!                rndy1=rndmin + (rndmax - rndmin)*rndx
!                call random_number(rndx)
!                rndy2=rndmin + (rndmax - rndmin)*rndx
!                !u*(m1)
!                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy1)
!                duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy1)
!                dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy1)
!                duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy1)
!                !u(m1) 
!                b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy2)
!                b_duz_dxi = b_duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy2)
!                b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy2)
!                b_duz_dgamma = b_duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy2)
!                !--------------------------------------------------------------------------
    
          enddo

          ! 1.2 lucas, to get xixl, xizl, gammaxl, gammazl  
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          ! 1.3 lucas, to get dux_dxl, dux_dzl, duz_dxl, duz_dzl, and also the back counterparts
          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

          iglob = ibool(i,j,ispec)
          ! 1.4 lucas dsxx, dsxz, dszz, the back counterparts, and kappa_k, mu_k
          ! isotropic kernel contributions
          if (P_SV) then
            ! P-SV waves
            dsxx =  dux_dxl
            dsxz = HALF * (duz_dxl + dux_dzl)
            dszz =  duz_dzl

            b_dsxx =  b_dux_dxl
            b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
            b_dszz =  b_duz_dzl

            kappa_k(iglob) = (dsxx + dszz) *  (b_dsxx + b_dszz)
            mu_k(iglob) = dsxx * b_dsxx + dszz * b_dszz + &
                          2._CUSTOM_REAL * dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob) ! lucas, here should be 4/9 not 1/3.

            !lucas added for approximate hessian--------------------------------------------------------- 
            !kappa_k(iglob) = (b_dsxx + b_dszz) *  (b_dsxx + b_dszz)
            !mu_k(iglob) = b_dsxx * b_dsxx + b_dszz * b_dszz + &
            !              2._CUSTOM_REAL * b_dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob)
            ! --------------------------------------------------------------------------------------------

          else
            ! SH (membrane) waves
            mu_k(iglob) = dux_dxl * b_dux_dxl + dux_dzl * b_dux_dzl

            ! lucas added for approximate hessian-----------------------
            ! mu_k(iglob) = b_dux_dxl * b_dux_dxl + b_dux_dzl * b_dux_dzl
            ! ----------------------------------------------------------

          endif

          ! Voigt kernels, e.g., see Sieminski, 2007a,b
          if (ispec_is_anisotropic(ispec)) then
            c11_k(iglob) = dux_dxl*b_dux_dxl
            c13_k(iglob) = dux_dxl*b_duz_dzl + duz_dzl*b_dux_dxl
            c15_k(iglob) = 2*(dux_dxl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_dux_dxl)
            c33_k(iglob) = duz_dzl*b_duz_dzl
            c35_k(iglob) = 2*(duz_dzl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_duz_dzl)
            c55_k(iglob) = 4*HALF*(dux_dzl+duz_dxl)*HALF*(b_dux_dzl+b_duz_dxl)
          endif
        enddo;enddo
      endif
    enddo

    do iglob = 1,nglob
      ! rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) + accel_elastic(2,iglob)*b_displ_elastic(2,iglob) ! lucas, where 1 means x, and 2 means z component.
       rho_k(iglob) =  displ_elastic(1,iglob)*b_accel_elastic(1,iglob) + displ_elastic(2,iglob)*b_accel_elastic(2,iglob) ! lucas, where 1 means x, and 2 means z component.

       !lucas added for approximate Hessian------------------------------------------------------------------------------ 
       !rho_k(iglob) =  b_veloc_elastic(1,iglob)*b_veloc_elastic(1,iglob) + b_veloc_elastic(2,iglob)*b_veloc_elastic(2,iglob)
       ! ----------------------------------------------------------------------------------------------------------------

       !lucas: for test----------------------------------
       iglob_source_lucas=ibool(2,2,ispec_selected_rec(1))
       if(myrank==islice_selected_rec(1)) then
        if(iglob==iglob_source_lucas) then
         if(mod(it,500)==0) then 
           write(6,*) 'in compute_kernels to print adjoint and forward fields:'
           write(6,*) 'ispec_selected_source(1) = ',ispec_selected_rec(1)
           write(6,*) 'iglob_source_lucas = ',iglob_source_lucas
           write(6,*) 'displ_elastic(2,iglob_source_lucas) = u*(m) = ', &
                       displ_elastic(2,iglob_source_lucas)  
           write(6,*) 'b_displ_elastic(2,iglob_source_lucas) = u(m) = ', & 
                       b_displ_elastic(2,iglob_source_lucas) 
         endif
        endif
       endif
       !lucas: for test--------------------------------------
    enddo

  else
    ! updates kernels on GPU
    call compute_kernels_elastic_cuda(Mesh_pointer,deltatf)
  endif

 !2 .lucas, to get kernels (rho_kl, mu_kl, kappa_kl) at each iglob
 do ispec = 1, nspec
    if (ispec_is_elastic(ispec)) then
      ! isotropic kernels
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)

          ! for parameterization (rho,mu,kappa): "primary" kernels
          ! density kernel
          rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rho_k(iglob) ! lucas, the rho_kl() initializes at prepare_timerun.f90, in particuler in prepare_timerun_adjoint().
          ! shear modulus kernel
          mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) -  mu_k(iglob)
          ! bulk modulus kernel
          kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) -  kappa_k(iglob)

        enddo
      enddo
      ! Voigt kernels, e.g., see Sieminski, 2007a,b
      if (ispec_is_anisotropic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rho_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c11_kl(i,j,ispec) = c11_kl(i,j,ispec) - c11_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c13_kl(i,j,ispec) = c13_kl(i,j,ispec) - c13_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c15_kl(i,j,ispec) = c15_kl(i,j,ispec) - c15_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c33_kl(i,j,ispec) = c33_kl(i,j,ispec) - c33_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c35_kl(i,j,ispec) = c35_kl(i,j,ispec) - c35_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c55_kl(i,j,ispec) = c55_kl(i,j,ispec) - c55_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + c11_kl(i,j,ispec) + &
                                 c13_kl(i,j,ispec) + c15_kl(i,j,ispec) + c33_kl(i,j,ispec) + &
                                 c35_kl(i,j,ispec) + c55_kl(i,j,ispec)
          enddo
        enddo
      endif


    endif
  enddo

  ! 3.lucas, multiply delta and parameeter at the last step since A0dt+A1dt+A2dt + ...+ Andt=(A0+A1+A2+...+An)dt, i.e., any constant at each time step can be done after the sum.
  !   lucas, we compute the kernels at a specfic it, which means that this if only be executed when the condition is satisfied.
  ! only at the last time step we multiply by delta and parameter value, it is not necessary to do it at each iteration

   if (NSTEP - it == mod(NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS)) then  !lucas:NSTEP_BETWEEN_COMPUTE_KERNELS=1 set in Par_file. mod(NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS) =0 
  !if (NSTEP-it == 2600) then ! lucas: 7400=10000-2600 ! lucas test, we save the result at 7400 step, not until the last step.
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        ! isotropic kernels
        do j = 1, NGLLZ
          do i = 1, NGLLX

            if (.not. assign_external_model) then
              rhol = density(1,kmato(ispec))
              mul = poroelastcoef(2,1,kmato(ispec))
              !if (AXISYM) then ! ABAB !!
              !Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
                kappal = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS * mul
              !else
              !  kappal = poroelastcoef(3,1,kmato(ispec)) - mul
              !endif
            else
              rhol = rhoext(i,j,ispec)
              mul = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
              !if (AXISYM) then ! ABAB !!
              ! Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
                kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - FOUR_THIRDS * mul
              !else
              !  kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - mul
              !endif
            endif

            ! for parameterization (rho,mu,kappa): "primary" kernels
            ! density kernel
            rho_kl(i,j,ispec) = rhol * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl(i,j,ispec) ! lucas see tromp et al. 2005
            ! shear modulus kernel
            mu_kl(i,j,ispec) =  TWO * mul * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl(i,j,ispec)
            ! bulk modulus kernel
            kappa_kl(i,j,ispec) = kappal * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl(i,j,ispec)

            ! for parameterization (rho,beta,alpha):
            ! rho prime kernel
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
            ! Vs kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
            beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - FOUR_THIRDS * mul/kappal * kappa_kl(i,j,ispec))
            ! Vp kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
            alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul/kappal) * kappa_kl(i,j,ispec)
            ! for bulk velocity c parameterization (rho,bulk_c,beta):
            bulk_c_kl(i,j,ispec) =  TWO * kappa_kl(i,j,ispec)
            bulk_beta_kl(i,j,ispec) =  TWO * mu_kl(i,j,ispec)

          enddo
        enddo
      endif ! elastic
    enddo !nspec loop
  !endif
  endif ! it == NSTEP

  end subroutine compute_kernels_el


!--------------------lucas, CTD-SEM--------------------------------------------------------------

  subroutine compute_kernels_el_Ha_Hb_Hc_Habc() ! where Hb=Hbm + Hbs(appro Hessian). included Frechet.
 !lucas, compute full Hessian hernels, Ha, Hb, Hc, and Habc=Ha+Hb+Hc, where each component of Hc is of each row of the Hc matrix.
 !lucas, the original variables/matrices used for Ha, and the variables/matrices with '_m2' used for Hbm, 
 !lucas, the variables/matrices with '_m1' used for Hbs means the approximate Hessian, 
 !lucas, Hc is computed based on Frechet and dm, that is to multiple the dm with the Frechet.
 !lucas, ref: Xie et al, Hessian Kernels on the fly by adjoint methods, GJI (2020).
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,HALF,TWO,FOUR_THIRDS

  use specfem_par, only: ispec_is_elastic,Full_Hessian_by_Wavefield_Stored, & !lucas, CTD-SEM
                         nglob,nspec,ibool,density,poroelastcoef,kmato,assign_external_model, & 
                         rhoext,vsext,vpext,rhoext_m2,vsext_m2,vpext_m2,deltat,P_SV, &
                         displ_elastic,b_displ_elastic,b_accel_elastic, & 
                         displ_elastic_m2,displ_elastic_m1,b_accel_elastic_m2,b_displ_elastic_m2, & !lucas
                         hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         GPU_MODE,it,NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS, &
                         no_backward_displ_buffer,no_backward_displ_buffer_fwd_du, &! lucas 
                         no_backward_accel_buffer_fwd_du,no_backward_displ_buffer_adj_du_s, & !lucas
                         no_backward_displ_buffer_adj_du_m, & !lucas
                         rho_k,mu_k,kappa_k,rho_k_Ha,mu_k_Ha,kappa_k_Ha,rho_k_Hbm,mu_k_Hbm,kappa_k_Hbm, & !lucas
                         rho_k_Hbs,mu_k_Hbs,kappa_k_Hbs,& !lucas
                         rho_kl,mu_kl,kappa_kl,rhop_kl,beta_kl,alpha_kl,& !bulk_c_kl,bulk_beta_kl, &
                         rho_kl_Ha,rho_kl_Hbm,rho_kl_Hbs,rho_kl_Hc,rho_kl_Habc,& !rho
                         mu_kl_Ha,mu_kl_Hbm,mu_kl_Hbs,mu_kl_Hc,mu_kl_Habc,& !mu
                         kappa_kl_Ha,kappa_kl_Hbm,kappa_kl_Hbs,kappa_kl_Hc,kappa_kl_Habc,& !kappa
                         rhop_kl_Ha,rhop_kl_Hbm,rhop_kl_Hbs,rhop_kl_Hc,rhop_kl_Habc,& !rhop
                         beta_kl_Ha,beta_kl_Hbm,beta_kl_Hbs,beta_kl_Hc,beta_kl_Habc,& !beta
                         alpha_kl_Ha,alpha_kl_Hbm,alpha_kl_Hbs,alpha_kl_Hc,alpha_kl_Habc !alpha
                         

  use specfem_par_gpu, only: Mesh_pointer,deltatf

  use specfem_par, only: c11_k,c13_k,c15_k,c33_k,c35_k,c55_k,ispec_is_anisotropic, &
                         rho_kl,c11_kl,c13_kl,c15_kl,c33_kl,c35_kl,c55_kl

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob
  !for Hc
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal
  real(kind=CUSTOM_REAL) :: rhol_m2,mul_m2,kappal_m2
  !for Ha
  real(kind=CUSTOM_REAL) :: dux_dxi_Ha,dux_dgamma_Ha,duz_dxi_Ha,duz_dgamma_Ha !i
  real(kind=CUSTOM_REAL) :: b_dux_dxi_Ha,b_dux_dgamma_Ha,b_duz_dxi_Ha,b_duz_dgamma_Ha 
  real(kind=CUSTOM_REAL) :: dux_dxl_Ha,dux_dzl_Ha,duz_dxl_Ha,duz_dzl_Ha  !l
  real(kind=CUSTOM_REAL) :: b_dux_dxl_Ha, b_dux_dzl_Ha,b_duz_dxl_Ha,b_duz_dzl_Ha  
  real(kind=CUSTOM_REAL) :: dsxx_Ha,dsxz_Ha,dszz_Ha
  real(kind=CUSTOM_REAL) :: b_dsxx_Ha,b_dsxz_Ha,b_dszz_Ha 
  
  !for Hbm
  real(kind=CUSTOM_REAL) :: dux_dxi_Hbm,dux_dgamma_Hbm,duz_dxi_Hbm,duz_dgamma_Hbm !i
  real(kind=CUSTOM_REAL) :: b_dux_dxi_Hbm,b_dux_dgamma_Hbm,b_duz_dxi_Hbm,b_duz_dgamma_Hbm 
  real(kind=CUSTOM_REAL) :: dux_dxl_Hbm,dux_dzl_Hbm,duz_dxl_Hbm,duz_dzl_Hbm  !l
  real(kind=CUSTOM_REAL) :: b_dux_dxl_Hbm,b_dux_dzl_Hbm,b_duz_dxl_Hbm,b_duz_dzl_Hbm 
  real(kind=CUSTOM_REAL) :: dsxx_Hbm,dsxz_Hbm,dszz_Hbm 
  real(kind=CUSTOM_REAL) :: b_dsxx_Hbm,b_dsxz_Hbm,b_dszz_Hbm 
  
  !for Hbs (for approximate Hessian)
  real(kind=CUSTOM_REAL) :: dux_dxi_Hbs,dux_dgamma_Hbs,duz_dxi_Hbs,duz_dgamma_Hbs !i
  real(kind=CUSTOM_REAL) :: b_dux_dxi_Hbs,b_dux_dgamma_Hbs,b_duz_dxi_Hbs,b_duz_dgamma_Hbs 
  real(kind=CUSTOM_REAL) :: dux_dxl_Hbs,dux_dzl_Hbs,duz_dxl_Hbs,duz_dzl_Hbs  !l
  real(kind=CUSTOM_REAL) :: b_dux_dxl_Hbs, b_dux_dzl_Hbs,b_duz_dxl_Hbs,b_duz_dzl_Hbs 
  real(kind=CUSTOM_REAL) :: dsxx_Hbs,dsxz_Hbs,dszz_Hbs 
  real(kind=CUSTOM_REAL) :: b_dsxx_Hbs,b_dsxz_Hbs,b_dszz_Hbs 
!  real(kind=CUSTOM_REAL) :: rndx,rndy1,rndy2,rndy3,rndy4,rndy5 !lucas
 
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl
!  double precision :: rndmin, rndmax      ! lucas
!  rndmin = -0.1 !0.1=10%, 0.01=1%, 0.001=0.1%, 0.0001=0.01%
!  rndmax =  0.1

  ! 1.lucas, to get kappa_k, mu_k, and rho_k at each iglob
  ! elastic kernels
  if (.not. GPU_MODE) then
    ! updates kernels on CPU
    do ispec = 1,nspec
      if (ispec_is_elastic(ispec)) then
        do j = 1,NGLLZ; do i = 1,NGLLX
          ! derivative along x and along z
          ! for Hc, the orignal version
          dux_dxi = 0._CUSTOM_REAL
          duz_dxi = 0._CUSTOM_REAL
          dux_dgamma = 0._CUSTOM_REAL
          duz_dgamma = 0._CUSTOM_REAL
          b_dux_dxi = 0._CUSTOM_REAL
          b_duz_dxi = 0._CUSTOM_REAL
          b_dux_dgamma = 0._CUSTOM_REAL
          b_duz_dgamma = 0._CUSTOM_REAL
          ! for Ha
          dux_dxi_Ha = 0._CUSTOM_REAL
          duz_dxi_Ha = 0._CUSTOM_REAL
          dux_dgamma_Ha = 0._CUSTOM_REAL
          duz_dgamma_Ha = 0._CUSTOM_REAL
          b_dux_dxi_Ha = 0._CUSTOM_REAL
          b_duz_dxi_Ha = 0._CUSTOM_REAL
          b_dux_dgamma_Ha = 0._CUSTOM_REAL
          b_duz_dgamma_Ha = 0._CUSTOM_REAL
          ! for Hbm
          dux_dxi_Hbm = 0._CUSTOM_REAL
          duz_dxi_Hbm = 0._CUSTOM_REAL
          dux_dgamma_Hbm = 0._CUSTOM_REAL
          duz_dgamma_Hbm = 0._CUSTOM_REAL
          b_dux_dxi_Hbm = 0._CUSTOM_REAL
          b_duz_dxi_Hbm = 0._CUSTOM_REAL
          b_dux_dgamma_Hbm = 0._CUSTOM_REAL
          b_duz_dgamma_Hbm = 0._CUSTOM_REAL
          ! for Hbs
          dux_dxi_Hbs = 0._CUSTOM_REAL
          duz_dxi_Hbs = 0._CUSTOM_REAL
          dux_dgamma_Hbs = 0._CUSTOM_REAL
          duz_dgamma_Hbs = 0._CUSTOM_REAL
          b_dux_dxi_Hbs = 0._CUSTOM_REAL
          b_duz_dxi_Hbs = 0._CUSTOM_REAL
          b_dux_dgamma_Hbs = 0._CUSTOM_REAL
          b_duz_dgamma_Hbs = 0._CUSTOM_REAL

          ! 1.1 lucas, to get dux_dxi, duz_dxi, dux_dgamma, duz_dgamma and also the back conterparts
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX 
             !for Ha==========================================
             !lucas, used for du, used in m1
             !lucas: here b_displ_elastic is changed to du = no_backward_displ_buffer_fwd_du (see read_forward_array.f90 in details)
             if(Full_Hessian_by_Wavefield_Stored) then!-----------------
             !u*(m1)
             dux_dxi_Ha = dux_dxi_Ha + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
             duz_dxi_Ha = duz_dxi_Ha + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
             dux_dgamma_Ha = dux_dgamma_Ha + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
             duz_dgamma_Ha = duz_dgamma_Ha + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
             !du
             b_dux_dxi_Ha = b_dux_dxi_Ha + no_backward_displ_buffer_fwd_du(1,ibool(k,j,ispec))*hprime_xx(i,k)
             b_duz_dxi_Ha = b_duz_dxi_Ha + no_backward_displ_buffer_fwd_du(2,ibool(k,j,ispec))*hprime_xx(i,k)
             b_dux_dgamma_Ha = b_dux_dgamma_Ha + no_backward_displ_buffer_fwd_du(1,ibool(i,k,ispec))*hprime_zz(j,k)
             b_duz_dgamma_Ha = b_duz_dgamma_Ha + no_backward_displ_buffer_fwd_du(2,ibool(i,k,ispec))*hprime_zz(j,k)
             else ! for CTD_SEM, without store the entire fields
             !u*(m1)
             dux_dxi_Ha = dux_dxi_Ha + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k) 
             duz_dxi_Ha = duz_dxi_Ha + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
             dux_dgamma_Ha = dux_dgamma_Ha + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
             duz_dgamma_Ha = duz_dgamma_Ha + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
             !du                                                                                        
             b_dux_dxi_Ha = b_dux_dxi_Ha + &
              (b_displ_elastic_m2(1,ibool(k,j,ispec))-b_displ_elastic(1,ibool(k,j,ispec)))*hprime_xx(i,k) ! du=u(m+dm)-d(m)
             b_duz_dxi_Ha = b_duz_dxi_Ha + &
              (b_displ_elastic_m2(2,ibool(k,j,ispec))-b_displ_elastic(2,ibool(k,j,ispec)))*hprime_xx(i,k)
             b_dux_dgamma_Ha = b_dux_dgamma_Ha + &
              (b_displ_elastic_m2(1,ibool(i,k,ispec))-b_displ_elastic(1,ibool(i,k,ispec)))*hprime_zz(j,k)
             b_duz_dgamma_Ha = b_duz_dgamma_Ha + &
              (b_displ_elastic_m2(2,ibool(i,k,ispec))-b_displ_elastic(2,ibool(i,k,ispec)))*hprime_zz(j,k)
             
!               !====================================================================================================
!               call random_number(rndx)
!               ! get a rondom number between [rndmin rndmax]
!               rndy1=rndmin + (rndmax - rndmin)*rndx
!               call random_number(rndx)
!               rndy2=rndmin + (rndmax - rndmin)*rndx
!               call random_number(rndx)
!               rndy3=rndmin + (rndmax - rndmin)*rndx
!               call random_number(rndx)
!               rndy4=rndmin + (rndmax - rndmin)*rndx
!               call random_number(rndx)
!               rndy5=rndmin + (rndmax - rndmin)*rndx

               !------------------------------------------Ha
               ! update u* and du for checking accurcy for mimicing field
               ! compression method
!               !u*
!               dux_dxi_Ha = dux_dxi_Ha + &
!               displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy1) 
!               duz_dxi_Ha = duz_dxi_Ha + &
!               displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy1)
!               dux_dgamma_Ha = dux_dgamma_Ha + &
!               displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy1)
!               duz_dgamma_Ha = duz_dgamma_Ha + &
!               displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy1) 
               !------
               !du  
!               b_dux_dxi_Ha = b_dux_dxi_Ha + &
!               (b_displ_elastic_m2(1,ibool(k,j,ispec))*(1+rndy2) - b_displ_elastic(1,ibool(k,j,ispec))*(1+rndy3))*hprime_xx(i,k)
!              ! du=u(m+dm)-d(m)
!               b_duz_dxi_Ha = b_duz_dxi_Ha + &
!               (b_displ_elastic_m2(2,ibool(k,j,ispec))*(1+rndy2) - b_displ_elastic(2,ibool(k,j,ispec))*(1+rndy3))*hprime_xx(i,k)
!               b_dux_dgamma_Ha = b_dux_dgamma_Ha + &
!               (b_displ_elastic_m2(1,ibool(i,k,ispec))*(1+rndy2) - b_displ_elastic(1,ibool(i,k,ispec))*(1+rndy3))*hprime_zz(j,k)
!               b_duz_dgamma_Ha = b_duz_dgamma_Ha + &
!               (b_displ_elastic_m2(2,ibool(i,k,ispec))*(1+rndy2) - b_displ_elastic(2,ibool(i,k,ispec))*(1+rndy3))*hprime_zz(j,k)               
               !-------------------------------------------------

            endif  ! Ha

            !for Hbm=================================================
            !lucas, used for du*=u*_{m}(m2)-u*(m1), adjoint source for both are adjsrc(m1)
            !lucas: here displ_elastic is changed to du* = no_backward_displ_buffer_adj_du_m (see read_forward_array.f90 in details)
            if(Full_Hessian_by_Wavefield_Stored) then
            !du*
            dux_dxi_Hbm = dux_dxi_Hbm + no_backward_displ_buffer_adj_du_m(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi_Hbm = duz_dxi_Hbm + no_backward_displ_buffer_adj_du_m(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma_Hbm = dux_dgamma_Hbm + no_backward_displ_buffer_adj_du_m(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma_Hbm = duz_dgamma_Hbm + no_backward_displ_buffer_adj_du_m(2,ibool(i,k,ispec))*hprime_zz(j,k)
            !u(m1)
            b_dux_dxi_Hbm = b_dux_dxi_Hbm + no_backward_displ_buffer(1,ibool(k,j,ispec))*hprime_xx(i,k) ! b_displ_elastic
            b_duz_dxi_Hbm = b_duz_dxi_Hbm + no_backward_displ_buffer(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma_Hbm = b_dux_dgamma_Hbm + no_backward_displ_buffer(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma_Hbm = b_duz_dgamma_Hbm + no_backward_displ_buffer(2,ibool(i,k,ispec))*hprime_zz(j,k)
            else ! lucas, CTD-SEM
            !du*
            dux_dxi_Hbm = dux_dxi_Hbm + &
             (displ_elastic_m2(1,ibool(k,j,ispec))-displ_elastic(1,ibool(k,j,ispec)))*hprime_xx(i,k) !du*=u*_{m}(m2)-u*(m1), both with adjsrc(m1)
            duz_dxi_Hbm = duz_dxi_Hbm + & 
             (displ_elastic_m2(2,ibool(k,j,ispec))-displ_elastic(2,ibool(k,j,ispec)))*hprime_xx(i,k)
            dux_dgamma_Hbm=dux_dgamma_Hbm + &
             (displ_elastic_m2(1,ibool(i,k,ispec))-displ_elastic(1,ibool(i,k,ispec)))*hprime_zz(j,k)
            duz_dgamma_Hbm=duz_dgamma_Hbm + &
             (displ_elastic_m2(2,ibool(i,k,ispec))-displ_elastic(2,ibool(i,k,ispec)))*hprime_zz(j,k)
            
            !u(m1)
            b_dux_dxi_Hbm = b_dux_dxi_Hbm + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k) !u(m1)
            b_duz_dxi_Hbm = b_duz_dxi_Hbm + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma_Hbm = b_dux_dgamma_Hbm + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma_Hbm = b_duz_dgamma_Hbm + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
 
               !---------------------------------------------------Hbm
!               ! update du* and u for check accuracy for mimicing field compression method 
!               !du*
!                dux_dxi_Hbm = dux_dxi_Hbm + &
!               (displ_elastic_m2(1,ibool(k,j,ispec))*(1+rndy4)-displ_elastic(1,ibool(k,j,ispec))*(1+rndy1))*hprime_xx(i,k)
!               !du*=u*_{m}(m2)-u*(m1), both with adjsrc(m1)
!               duz_dxi_Hbm = duz_dxi_Hbm + & 
!               (displ_elastic_m2(2,ibool(k,j,ispec))*(1+rndy4)-displ_elastic(2,ibool(k,j,ispec))*(1+rndy1))*hprime_xx(i,k)
!               dux_dgamma_Hbm=dux_dgamma_Hbm + &
!               (displ_elastic_m2(1,ibool(i,k,ispec))*(1+rndy4)-displ_elastic(1,ibool(i,k,ispec))*(1+rndy1))*hprime_zz(j,k)
!               duz_dgamma_Hbm=duz_dgamma_Hbm + &
!               (displ_elastic_m2(2,ibool(i,k,ispec))*(1+rndy4)-displ_elastic(2,ibool(i,k,ispec))*(1+rndy1))*hprime_zz(j,k)
!               !u
!               b_dux_dxi_Hbm = b_dux_dxi_Hbm + &
!                               b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy3)
!                              !u(m1)
!               b_duz_dxi_Hbm = b_duz_dxi_Hbm + &
!                               b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy3)
!               b_dux_dgamma_Hbm = b_dux_dgamma_Hbm + &
!                               b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy3)
!               b_duz_dgamma_Hbm = b_duz_dgamma_Hbm + &
!                               b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy3)
!               !-------------------------------------------------------

            endif !Hbm

            !for Hbs=================================================
            !lucas, used for d_{s}u*=u*_{s}(m1)-u*(m1), where adjsrc(m2) for u*_{s}(m1), and adjsrc(m1) for u*(m1)
            if(Full_Hessian_by_Wavefield_Stored) then ! lucas, need to do here
            !du*
            dux_dxi_Hbs = dux_dxi_Hbs + no_backward_displ_buffer_adj_du_s(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi_Hbs = duz_dxi_Hbs + no_backward_displ_buffer_adj_du_s(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma_Hbs = dux_dgamma_Hbs + no_backward_displ_buffer_adj_du_s(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma_Hbs = duz_dgamma_Hbs + no_backward_displ_buffer_adj_du_s(2,ibool(i,k,ispec))*hprime_zz(j,k)
            !u(m1)
            b_dux_dxi_Hbs = b_dux_dxi_Hbs + no_backward_displ_buffer(1,ibool(k,j,ispec))*hprime_xx(i,k) ! b_displ_elastic
            b_duz_dxi_Hbs = b_duz_dxi_Hbs + no_backward_displ_buffer(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma_Hbs = b_dux_dgamma_Hbs + no_backward_displ_buffer(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma_Hbs = b_duz_dgamma_Hbs + no_backward_displ_buffer(2,ibool(i,k,ispec))*hprime_zz(j,k)
            else ! lucas, CTD-SEM
            !du*
            dux_dxi_Hbs = dux_dxi_Hbs + &
             (displ_elastic_m1(1,ibool(k,j,ispec))-displ_elastic(1,ibool(k,j,ispec)))*hprime_xx(i,k) !du*=u*_{s}(m1)-u*(m1), where adjsrc(m2) for u*_{s}(m1)
            duz_dxi_Hbs = duz_dxi_Hbs + &
             (displ_elastic_m1(2,ibool(k,j,ispec))-displ_elastic(2,ibool(k,j,ispec)))*hprime_xx(i,k)
            dux_dgamma_Hbs=dux_dgamma_Hbs + &
             (displ_elastic_m1(1,ibool(i,k,ispec))-displ_elastic(1,ibool(i,k,ispec)))*hprime_zz(j,k)
            duz_dgamma_Hbs=duz_dgamma_Hbs + &
             (displ_elastic_m1(2,ibool(i,k,ispec))-displ_elastic(2,ibool(i,k,ispec)))*hprime_zz(j,k)
            !u(m1)
            b_dux_dxi_Hbs = b_dux_dxi_Hbs + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k) !u(m1)
            b_duz_dxi_Hbs = b_duz_dxi_Hbs + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma_Hbs = b_dux_dgamma_Hbs + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma_Hbs = b_duz_dgamma_Hbs + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k) 

!               !--------------------------------------------Hbs
!               ! update du* and u for check accuracy for mimicing field
!               ! compression method 
!               !du*
!                dux_dxi_Hbs = dux_dxi_Hbs + &
!               (displ_elastic_m1(1,ibool(k,j,ispec))*(1+rndy5)-displ_elastic(1,ibool(k,j,ispec))*(1+rndy1))*hprime_xx(i,k)
!               !du*=u*_{s}(m2)-u*(m1)
!               duz_dxi_Hbs = duz_dxi_Hbs + &
!               (displ_elastic_m1(2,ibool(k,j,ispec))*(1+rndy5)-displ_elastic(2,ibool(k,j,ispec))*(1+rndy1))*hprime_xx(i,k)
!               dux_dgamma_Hbs=dux_dgamma_Hbs + &
!               (displ_elastic_m1(1,ibool(i,k,ispec))*(1+rndy5)-displ_elastic(1,ibool(i,k,ispec))*(1+rndy1))*hprime_zz(j,k)
!               duz_dgamma_Hbs=duz_dgamma_Hbs + &
!               (displ_elastic_m1(2,ibool(i,k,ispec))*(1+rndy5)-displ_elastic(2,ibool(i,k,ispec))*(1+rndy1))*hprime_zz(j,k)
!               !u
!               b_dux_dxi_Hbs = b_dux_dxi_Hbs + &
!                               b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy3)
!                              !u(m1)
!               b_duz_dxi_Hbs = b_duz_dxi_Hbs + &
!                               b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy3)
!               b_dux_dgamma_Hbs = b_dux_dgamma_Hbs + &
!                               b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy3)
!               b_duz_dgamma_Hbs = b_duz_dgamma_Hbs + &
!                               b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy3)
!               !------------------------------------------------ 

            endif !Hbs

            ! for Hc, the original version========================================
            if(Full_Hessian_by_Wavefield_Stored) then
            !u*(m1)
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
            !u(m1) 
            b_dux_dxi = b_dux_dxi + no_backward_displ_buffer(1,ibool(k,j,ispec))*hprime_xx(i,k)
            b_duz_dxi = b_duz_dxi + no_backward_displ_buffer(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma = b_dux_dgamma + no_backward_displ_buffer(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma = b_duz_dgamma + no_backward_displ_buffer(2,ibool(i,k,ispec))*hprime_zz(j,k)
            else !lucas, CTD-SEM
            !u*(m1)
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
            !u(m1) 
            b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            b_duz_dxi = b_duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma = b_duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

!              !----------------------------------------------Hc
!               ! update u* and u for check accuracy for mimicing field
!               !u*
!               dux_dxi = dux_dxi + &
!               displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy1)
!               duz_dxi = duz_dxi + &
!               displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy1)
!               dux_dgamma = dux_dgamma + &
!               displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy1)
!               duz_dgamma = duz_dgamma + &
!               displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy1)
!               !u(m1) 
!               b_dux_dxi = b_dux_dxi + &
!               b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy3)
!               b_duz_dxi = b_duz_dxi + &
!               b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)*(1+rndy3)
!               b_dux_dgamma = b_dux_dgamma + &
!               b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy3)
!               b_duz_dgamma = b_duz_dgamma + &
!               b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)*(1+rndy3)            
!              !------------------------------------------------

            endif !Hc

!           !=================================================================================================== 
           
          enddo ! end k for GLL

          ! 1.2 lucas, to get xixl, xizl, gammaxl, gammazl  
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          ! 1.3 lucas, to get dux_dxl, dux_dzl, duz_dxl, duz_dzl, and also the back counterparts
          ! derivatives of displacement
          ! for Hc
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl
          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
          ! for Ha
          dux_dxl_Ha = dux_dxi_Ha*xixl + dux_dgamma_Ha*gammaxl
          dux_dzl_Ha = dux_dxi_Ha*xizl + dux_dgamma_Ha*gammazl
          duz_dxl_Ha = duz_dxi_Ha*xixl + duz_dgamma_Ha*gammaxl
          duz_dzl_Ha = duz_dxi_Ha*xizl + duz_dgamma_Ha*gammazl

          b_dux_dxl_Ha = b_dux_dxi_Ha*xixl + b_dux_dgamma_Ha*gammaxl
          b_dux_dzl_Ha = b_dux_dxi_Ha*xizl + b_dux_dgamma_Ha*gammazl
          b_duz_dxl_Ha = b_duz_dxi_Ha*xixl + b_duz_dgamma_Ha*gammaxl
          b_duz_dzl_Ha = b_duz_dxi_Ha*xizl + b_duz_dgamma_Ha*gammazl
          !for Hbm
          dux_dxl_Hbm = dux_dxi_Hbm*xixl + dux_dgamma_Hbm*gammaxl
          dux_dzl_Hbm = dux_dxi_Hbm*xizl + dux_dgamma_Hbm*gammazl
          duz_dxl_Hbm = duz_dxi_Hbm*xixl + duz_dgamma_Hbm*gammaxl
          duz_dzl_Hbm = duz_dxi_Hbm*xizl + duz_dgamma_Hbm*gammazl

          b_dux_dxl_Hbm = b_dux_dxi_Hbm*xixl + b_dux_dgamma_Hbm*gammaxl
          b_dux_dzl_Hbm = b_dux_dxi_Hbm*xizl + b_dux_dgamma_Hbm*gammazl
          b_duz_dxl_Hbm = b_duz_dxi_Hbm*xixl + b_duz_dgamma_Hbm*gammaxl
          b_duz_dzl_Hbm = b_duz_dxi_Hbm*xizl + b_duz_dgamma_Hbm*gammazl

          !for Hbs
          dux_dxl_Hbs = dux_dxi_Hbs*xixl + dux_dgamma_Hbs*gammaxl
          dux_dzl_Hbs = dux_dxi_Hbs*xizl + dux_dgamma_Hbs*gammazl
          duz_dxl_Hbs = duz_dxi_Hbs*xixl + duz_dgamma_Hbs*gammaxl
          duz_dzl_Hbs = duz_dxi_Hbs*xizl + duz_dgamma_Hbs*gammazl

          b_dux_dxl_Hbs = b_dux_dxi_Hbs*xixl + b_dux_dgamma_Hbs*gammaxl
          b_dux_dzl_Hbs = b_dux_dxi_Hbs*xizl + b_dux_dgamma_Hbs*gammazl
          b_duz_dxl_Hbs = b_duz_dxi_Hbs*xixl + b_duz_dgamma_Hbs*gammaxl
          b_duz_dzl_Hbs = b_duz_dxi_Hbs*xizl + b_duz_dgamma_Hbs*gammazl

          iglob = ibool(i,j,ispec)
          ! 1.4 lucas dsxx, dsxz, dszz, the back counterparts, and kappa_k, mu_k
          ! isotropic kernel contributions
          if (P_SV) then
            ! P-SV waves
            !for Hc
            dsxx =  dux_dxl
            dsxz = HALF * (duz_dxl + dux_dzl)
            dszz =  duz_dzl

            b_dsxx =  b_dux_dxl
            b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
            b_dszz =  b_duz_dzl

           !lucas add initilize
            mu_k(iglob)=0._CUSTOM_REAL
            kappa_k(iglob)=0._CUSTOM_REAL
            mu_k_Ha(iglob)=0._CUSTOM_REAL
            kappa_k_Ha(iglob)=0._CUSTOM_REAL
            mu_k_Hbm(iglob)=0._CUSTOM_REAL
            kappa_k_Hbm(iglob)=0._CUSTOM_REAL
            mu_k_Hbs(iglob)=0._CUSTOM_REAL
            kappa_k_Hbs(iglob)=0._CUSTOM_REAL


            kappa_k(iglob) = (dsxx + dszz) *  (b_dsxx + b_dszz)
            mu_k(iglob) = dsxx * b_dsxx + dszz * b_dszz + &
                          2._CUSTOM_REAL * dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob)
            !for Ha
            dsxx_Ha =  dux_dxl_Ha
            dsxz_Ha = HALF * (duz_dxl_Ha + dux_dzl_Ha)
            dszz_Ha =  duz_dzl_Ha

            b_dsxx_Ha =  b_dux_dxl_Ha
            b_dsxz_Ha = HALF * (b_duz_dxl_Ha + b_dux_dzl_Ha)
            b_dszz_Ha =  b_duz_dzl_Ha

            kappa_k_Ha(iglob) = (dsxx_Ha + dszz_Ha) *  (b_dsxx_Ha + b_dszz_Ha)
            mu_k_Ha(iglob) = dsxx_Ha * b_dsxx_Ha + dszz_Ha * b_dszz_Ha + &
                          2._CUSTOM_REAL * dsxz_Ha * b_dsxz_Ha - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k_Ha(iglob)
            !for Hbm
            dsxx_Hbm =  dux_dxl_Hbm
            dsxz_Hbm = HALF * (duz_dxl_Hbm + dux_dzl_Hbm)
            dszz_Hbm =  duz_dzl_Hbm

            b_dsxx_Hbm =  b_dux_dxl_Hbm
            b_dsxz_Hbm = HALF * (b_duz_dxl_Hbm + b_dux_dzl_Hbm)
            b_dszz_Hbm =  b_duz_dzl_Hbm

            kappa_k_Hbm(iglob) = (dsxx_Hbm + dszz_Hbm) *  (b_dsxx_Hbm + b_dszz_Hbm)
            mu_k_Hbm(iglob) = dsxx_Hbm * b_dsxx_Hbm + dszz_Hbm * b_dszz_Hbm + &
                          2._CUSTOM_REAL * dsxz_Hbm * b_dsxz_Hbm - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k_Hbm(iglob)

            !for Hbs
            dsxx_Hbs =  dux_dxl_Hbs
            dsxz_Hbs = HALF * (duz_dxl_Hbs + dux_dzl_Hbs)
            dszz_Hbs =  duz_dzl_Hbs

            b_dsxx_Hbs =  b_dux_dxl_Hbs
            b_dsxz_Hbs = HALF * (b_duz_dxl_Hbs + b_dux_dzl_Hbs)
            b_dszz_Hbs =  b_duz_dzl_Hbs

            kappa_k_Hbs(iglob) = (dsxx_Hbs + dszz_Hbs) *  (b_dsxx_Hbs + b_dszz_Hbs)
            mu_k_Hbs(iglob) = dsxx_Hbs * b_dsxx_Hbs + dszz_Hbs * b_dszz_Hbs + &
                          2._CUSTOM_REAL * dsxz_Hbs * b_dsxz_Hbs - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k_Hbs(iglob)

          else
            ! SH (membrane) waves
            ! for Hc
            mu_k(iglob) = dux_dxl * b_dux_dxl + dux_dzl * b_dux_dzl
            ! for Ha
            mu_k_Ha(iglob) = dux_dxl_Ha * b_dux_dxl_Ha + dux_dzl_Ha * b_dux_dzl_Ha
            ! for Hbm
            mu_k_Hbm(iglob) = dux_dxl_Hbm * b_dux_dxl_Hbm + dux_dzl_Hbm * b_dux_dzl_Hbm
            ! for Hbs
            mu_k_Hbs(iglob) = dux_dxl_Hbs * b_dux_dxl_Hbs + dux_dzl_Hbs * b_dux_dzl_Hbs
          endif

          ! Voigt kernels, e.g., see Sieminski, 2007a,b
          if (ispec_is_anisotropic(ispec)) then
            c11_k(iglob) = dux_dxl*b_dux_dxl
            c13_k(iglob) = dux_dxl*b_duz_dzl + duz_dzl*b_dux_dxl
            c15_k(iglob) = 2*(dux_dxl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_dux_dxl)
            c33_k(iglob) = duz_dzl*b_duz_dzl
            c35_k(iglob) = 2*(duz_dzl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_duz_dzl)
            c55_k(iglob) = 4*HALF*(dux_dzl+duz_dxl)*HALF*(b_dux_dzl+b_duz_dxl)
          endif
        enddo;enddo ! end i, j
      endif ! end, elastic 
    enddo ! end nspec

    do iglob = 1,nglob ! lucas, here Hc, Ha, Hb are used in stored method.
       !lucas add initialize
       rho_k(iglob)=0._CUSTOM_REAL
       rho_k_Ha(iglob)=0._CUSTOM_REAL
       rho_k_Hbm(iglob)=0._CUSTOM_REAL
       rho_k_Hbs(iglob)=0._CUSTOM_REAL
       !lucas, the accelaration fields need to be updated here for density kernel caculations in case for imitating the
       !compression method by adding random noise. usualy do not need to do this
       !since only the K_alpha or K_beta are shown in publications, BBB

       ! for Hc, the same for stored-method and CTD-SEM
       !rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) + accel_elastic(2,iglob)*b_displ_elastic(2,iglob) ! lucas, where 1 means x, and 2 means z component.
       rho_k(iglob) =     displ_elastic(1,iglob)*b_accel_elastic(1,iglob) + displ_elastic(2,iglob)*b_accel_elastic(2,iglob) ! lucas, b_accel_elastic = read from disk or by b_simulation
       if(Full_Hessian_by_Wavefield_Stored) then!-----------
       ! for Ha (du,u*)
       rho_k_Ha(iglob) =  displ_elastic(1,iglob)*no_backward_accel_buffer_fwd_du(1,iglob) + &  !lucas, u*, du_accel
                          displ_elastic(2,iglob)*no_backward_accel_buffer_fwd_du(2,iglob)
       ! for Hbm (u,du*)
       rho_k_Hbm(iglob) =  no_backward_displ_buffer_adj_du_m(1,iglob)*b_accel_elastic(1,iglob) + & !lucas, du*, u_accel
                          no_backward_displ_buffer_adj_du_m(2,iglob)*b_accel_elastic(2,iglob)  
       ! for Hbs (u,du*)
       rho_k_Hbs(iglob) =  no_backward_displ_buffer_adj_du_s(1,iglob)*b_accel_elastic(1,iglob) + & !lucas, du*, u_accel
                          no_backward_displ_buffer_adj_du_s(2,iglob)*b_accel_elastic(2,iglob)  
       else !lucas, CTD-SEM
       ! for Ha (du,u*)
       rho_k_Ha(iglob) =  displ_elastic(1,iglob)*(b_accel_elastic_m2(1,iglob)-b_accel_elastic(1,iglob)) + &  !lucas, u*, du_accel
                          displ_elastic(2,iglob)*(b_accel_elastic_m2(2,iglob)-b_accel_elastic(2,iglob))
       ! for Hbm (u,du*)
       rho_k_Hbm(iglob) = (displ_elastic_m2(1,iglob)-displ_elastic(1,iglob))*b_accel_elastic(1,iglob) + & !lucas, du*_{m}, u_accel
                          (displ_elastic_m2(2,iglob)-displ_elastic(2,iglob))*b_accel_elastic(2,iglob)  

       ! for Hbs (u,du*)
       rho_k_Hbs(iglob) = (displ_elastic_m1(1,iglob)-displ_elastic(1,iglob))*b_accel_elastic(1,iglob) + & !lucas, du*_{s}, u_accel
                          (displ_elastic_m1(2,iglob)-displ_elastic(2,iglob))*b_accel_elastic(2,iglob)  
       endif!---------
    enddo

  else
    ! updates kernels on GPU
    call compute_kernels_elastic_cuda(Mesh_pointer,deltatf)
  endif

 !2 .lucas, to get kernels (rho_kl, mu_kl, kappa_kl) at each iglob
 do ispec = 1, nspec
    if (ispec_is_elastic(ispec)) then
      ! isotropic kernels
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)

          ! for parameterization (rho,mu,kappa): "primary" kernels
          !for Hc, it is composed of frechet
          ! density kernel
          rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rho_k(iglob) 
          ! shear modulus kernel
          mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) -  mu_k(iglob)
          ! bulk modulus kernel
          kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) -  kappa_k(iglob)
          !for Ha-----------
          ! density kernel
          rho_kl_Ha(i,j,ispec) = rho_kl_Ha(i,j,ispec) - rho_k_Ha(iglob) 
          ! shear modulus kernel
          mu_kl_Ha(i,j,ispec) =  mu_kl_Ha(i,j,ispec) -  mu_k_Ha(iglob)
          ! bulk modulus kernel
          kappa_kl_Ha(i,j,ispec) = kappa_kl_Ha(i,j,ispec) -  kappa_k_Ha(iglob)
          !for Hbm-----------
          ! density kernel
          rho_kl_Hbm(i,j,ispec) = rho_kl_Hbm(i,j,ispec) - rho_k_Hbm(iglob) 
          ! shear modulus kernel
          mu_kl_Hbm(i,j,ispec) =  mu_kl_Hbm(i,j,ispec) -  mu_k_Hbm(iglob)
          ! bulk modulus kernel
          kappa_kl_Hbm(i,j,ispec) = kappa_kl_Hbm(i,j,ispec) -  kappa_k_Hbm(iglob)

          !for Hbs-----------
          ! density kernel
          rho_kl_Hbs(i,j,ispec) = rho_kl_Hbs(i,j,ispec) - rho_k_Hbs(iglob) 
          ! shear modulus kernel
          mu_kl_Hbs(i,j,ispec) =  mu_kl_Hbs(i,j,ispec) -  mu_k_Hbs(iglob)
          ! bulk modulus kernel
          kappa_kl_Hbs(i,j,ispec) = kappa_kl_Hbs(i,j,ispec) -  kappa_k_Hbs(iglob)

        enddo
      enddo
      ! Voigt kernels, e.g., see Sieminski, 2007a,b
      if (ispec_is_anisotropic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rho_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c11_kl(i,j,ispec) = c11_kl(i,j,ispec) - c11_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c13_kl(i,j,ispec) = c13_kl(i,j,ispec) - c13_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c15_kl(i,j,ispec) = c15_kl(i,j,ispec) - c15_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c33_kl(i,j,ispec) = c33_kl(i,j,ispec) - c33_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c35_kl(i,j,ispec) = c35_kl(i,j,ispec) - c35_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c55_kl(i,j,ispec) = c55_kl(i,j,ispec) - c55_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + c11_kl(i,j,ispec) + &
                                 c13_kl(i,j,ispec) + c15_kl(i,j,ispec) + c33_kl(i,j,ispec) + &
                                 c35_kl(i,j,ispec) + c55_kl(i,j,ispec)
          enddo
        enddo
      endif


    endif
  enddo

  ! 3.lucas, multiply delta and parameeter at the last step since A0dt+A1dt+A2dt + ...+ Andt=(A0+A1+A2+...+An)dt, i.e., any constant at each time step can be done after the sum.
  !   lucas, we compute the kernels at a specfic it, which means tha this if only be execute when the condition is satisfied.
  ! only at the last time step we multiply by delta and parameter value, it is not necessary to do it at each iteration

   if (NSTEP - it == mod(NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS)) then  !lucas:NSTEP_BETWEEN_COMPUTE_KERNELS=1 set in Par_file. mod(NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS) =0
  ! lucas test, we save the result at 7400 step, not until the last step. 
  !if (NSTEP-it == 2600) then ! lucas: 7400=10000-2600
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        ! isotropic kernels
        do j = 1, NGLLZ
          do i = 1, NGLLX

            if (.not. assign_external_model) then
              rhol = density(1,kmato(ispec))
              mul = poroelastcoef(2,1,kmato(ispec))
              !if (AXISYM) then ! ABAB !!
              !Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
              kappal = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS * mul
              !else
              !kappal = poroelastcoef(3,1,kmato(ispec)) - mul
              !endif
            else
              ! for Ha and Hb, used the same model, i.e., m1
              rhol = rhoext(i,j,ispec)
              mul = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
              !if (AXISYM) then ! ABAB !!
              ! Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
              kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - FOUR_THIRDS * mul
              !else
              !kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - mul
              !endif
              !for Hc, load m1 and m2, where m1 is provied above, and m2=m1+dm is shown below.
              !for m2
              rhol_m2 = rhoext_m2(i,j,ispec)
              mul_m2 = rhoext_m2(i,j,ispec)*vsext_m2(i,j,ispec)*vsext_m2(i,j,ispec)
              kappal_m2 = rhoext_m2(i,j,ispec)*vpext_m2(i,j,ispec)*vpext_m2(i,j,ispec) - FOUR_THIRDS * mul_m2

            endif

            ! for parameterization (rho,mu,kappa): "primary" kernels
            !for Hc,-----------------------------------------------use m1 and m2---------
            ! density kernel
            rho_kl(i,j,ispec) = rhol * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl(i,j,ispec) 
            ! shear modulus kernel
            mu_kl(i,j,ispec) =  TWO * mul * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl(i,j,ispec)
            ! bulk modulus kernel
            kappa_kl(i,j,ispec) = kappal * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl(i,j,ispec)
            ! all componenets of Hc=0 when model is based on rho, mu, kappa.
            rho_kl_Hc(i,j,ispec) = 0._CUSTOM_REAL
            mu_kl_Hc(i,j,ispec) = 0._CUSTOM_REAL
            kappa_kl_Hc(i,j,ispec) = 0._CUSTOM_REAL

            ! for parameterization (rho,beta,alpha):
            ! rho prime kernel
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
            ! Vs kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
            beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - FOUR_THIRDS * mul/kappal * kappa_kl(i,j,ispec))
            ! Vp kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
            alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul/kappal) * kappa_kl(i,j,ispec)
            
            ! each items of the three below indicate each row of the Hc, where alpha_kl_Hc indicates the second row (Xie et al, GJI,2020)
            ! lucas,(rhoext, vpext, vsext) and (rhoext_m2, vpext_m2, vsext_m2) should > 0. 
            if(rhoext(i,j,ispec)>0 .and. vpext(i,j,ispec)>0 .and. vsext(i,j,ispec)>0) then ! lucas. only used for reading tomography from external files
            rhop_kl_Hc(i,j,ispec) = (1/rhoext(i,j,ispec))*alpha_kl(i,j,ispec)*(vpext_m2(i,j,ispec) - vpext(i,j,ispec)) + &
                                    (1/rhoext(i,j,ispec))*beta_kl(i,j,ispec)*(vsext_m2(i,j,ispec) - vsext(i,j,ispec))
            beta_kl_Hc(i,j,ispec) = (1/rhoext(i,j,ispec))*alpha_kl(i,j,ispec)*(rhoext_m2(i,j,ispec) - rhoext(i,j,ispec)) + &
                                    (1/vpext(i,j,ispec))*alpha_kl(i,j,ispec)*(vpext_m2(i,j,ispec)-vpext(i,j,ispec))
            alpha_kl_Hc(i,j,ispec) =(1/rhoext(i,j,ispec))*beta_kl(i,j,ispec)*(rhoext_m2(i,j,ispec) - rhoext(i,j,ispec)) + &
                                    (1/vsext(i,j,ispec))*beta_kl(i,j,ispec)*(vsext_m2(i,j,ispec) - vsext(i,j,ispec))
            else
            call stop_the_code('error: Vp Vs and Rho should be large than zero, i.e. for elastic only in this implementation,lucas')
            endif
            
            ! for bulk velocity c parameterization (rho,bulk_c,beta):
            ! bulk_c_kl(i,j,ispec) =  TWO * kappa_kl(i,j,ispec)
            ! bulk_beta_kl(i,j,ispec) =  TWO * mu_kl(i,j,ispec)


           ! for Ha--------------------------------------------------------------- use m1----------------------------
            ! density kernel
            rho_kl_Ha(i,j,ispec) = rhol * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl_Ha(i,j,ispec) 
            ! shear modulus kernel
            mu_kl_Ha(i,j,ispec) =  TWO * mul * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl_Ha(i,j,ispec)
            ! bulk modulus kernel
            kappa_kl_Ha(i,j,ispec) = kappal * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl_Ha(i,j,ispec)

            ! for parameterization (rho,beta,alpha):
            ! rho prime kernel
            rhop_kl_Ha(i,j,ispec) = rho_kl_Ha(i,j,ispec) + kappa_kl_Ha(i,j,ispec) + mu_kl_Ha(i,j,ispec)
            ! Vs kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
            beta_kl_Ha(i,j,ispec) = TWO * (mu_kl_Ha(i,j,ispec) - FOUR_THIRDS * mul/kappal * kappa_kl_Ha(i,j,ispec))
            ! Vp kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
            alpha_kl_Ha(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul/kappal) * kappa_kl_Ha(i,j,ispec)

           ! for Hbm-----------------------------------------------------------------use m1-------------------------------
            ! density kernel
            rho_kl_Hbm(i,j,ispec) = rhol * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl_Hbm(i,j,ispec) 
            ! shear modulus kernel
            mu_kl_Hbm(i,j,ispec) =  TWO * mul * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl_Hbm(i,j,ispec)
            ! bulk modulus kernel
            kappa_kl_Hbm(i,j,ispec) = kappal * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl_Hbm(i,j,ispec)

            ! for parameterization (rho,beta,alpha):
            ! rho prime kernel
            rhop_kl_Hbm(i,j,ispec) = rho_kl_Hbm(i,j,ispec) + kappa_kl_Hbm(i,j,ispec) + mu_kl_Hbm(i,j,ispec)
            ! Vs kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
            beta_kl_Hbm(i,j,ispec) = TWO * (mu_kl_Hbm(i,j,ispec) - FOUR_THIRDS * mul/kappal * kappa_kl_Hbm(i,j,ispec))
            ! Vp kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
            alpha_kl_Hbm(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul/kappal) * kappa_kl_Hbm(i,j,ispec)

            ! for Hbs-----------------------------------------------------------------use m1-------------------------------
            ! density kernel
            rho_kl_Hbs(i,j,ispec) = rhol * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl_Hbs(i,j,ispec) 
            ! shear modulus kernel
            mu_kl_Hbs(i,j,ispec) =  TWO * mul * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl_Hbs(i,j,ispec)
            ! bulk modulus kernel
            kappa_kl_Hbs(i,j,ispec) = kappal * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl_Hbs(i,j,ispec)

            ! for parameterization (rho,beta,alpha):
            ! rho prime kernel
            rhop_kl_Hbs(i,j,ispec) = rho_kl_Hbs(i,j,ispec) + kappa_kl_Hbs(i,j,ispec) + mu_kl_Hbs(i,j,ispec)
            ! Vs kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
            beta_kl_Hbs(i,j,ispec) = TWO * (mu_kl_Hbs(i,j,ispec) - FOUR_THIRDS * mul/kappal * kappa_kl_Hbs(i,j,ispec))
            ! Vp kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
            alpha_kl_Hbs(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul/kappal) * kappa_kl_Hbs(i,j,ispec)

            !for Habc=Ha+Hb+Hc, where Hb=Hbm+Hbs-----------------------------------------------------------------------------------------
         rho_kl_Habc(i,j,ispec)=rho_kl_Ha(i,j,ispec)+rho_kl_Hbm(i,j,ispec)+rho_kl_Hbs(i,j,ispec)+rho_kl_Hc(i,j,ispec)
         mu_kl_Habc(i,j,ispec)=mu_kl_Ha(i,j,ispec)+mu_kl_Hbm(i,j,ispec)+mu_kl_Hbs(i,j,ispec)+mu_kl_Hc(i,j,ispec)
         kappa_kl_Habc(i,j,ispec)=kappa_kl_Ha(i,j,ispec)+kappa_kl_Hbm(i,j,ispec)+kappa_kl_Hbs(i,j,ispec)+kappa_kl_Hc(i,j,ispec)

         rhop_kl_Habc(i,j,ispec)=rhop_kl_Ha(i,j,ispec)+rhop_kl_Hbm(i,j,ispec)+rhop_kl_Hbs(i,j,ispec)+rhop_kl_Hc(i,j,ispec)
         beta_kl_Habc(i,j,ispec)=beta_kl_Ha(i,j,ispec)+beta_kl_Hbm(i,j,ispec)+beta_kl_Hbs(i,j,ispec)+beta_kl_Hc(i,j,ispec)
         alpha_kl_Habc(i,j,ispec)=alpha_kl_Ha(i,j,ispec)+alpha_kl_Hbm(i,j,ispec)+alpha_kl_Hbs(i,j,ispec)+alpha_kl_Hc(i,j,ispec)

          enddo
        enddo
      endif ! elastic
    enddo !nspec loop
  endif ! it == NSTEP

  end subroutine compute_kernels_el_Ha_Hb_Hc_Habc





!--------------------compute_kernels_el_m2() for test only-------------------------------------------



!  subroutine compute_kernels_el_m2() !lucas, for simulation_type=3, the displ_elastic is obtained when the source is at the receiver, and b_displ_elastic read from last frame of the !forward
                                  !lucas, the first it here is the last step of the forward.
! elastic kernel calculations
! see e.g. Tromp et al. (2005)

!  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,HALF,TWO,FOUR_THIRDS

!  use specfem_par, only: ispec_is_elastic,rhop_kl, rho_k, & ! AXISYM,
!                         rho_kl_m2,mu_kl_m2,kappa_kl_m2,rhop_kl_m2,beta_kl_m2,alpha_kl_m2, & !lucas, CTD-SEM
!                         nglob,nspec,ibool, & !accel_elastic, & ! lucas added b_veloc_elastic for appro. hessian.  
!                         b_displ_elastic_m2, &  !lucas, CTD-SEM  
!                         density,poroelastcoef,kmato,assign_external_model, &
!                         rhoext_m2,vsext_m2,vpext_m2, & !lucas, CTD-SEM
!                         deltat,P_SV, &
!                         displ_elastic_m2, & !lucas, CTD-SEM
!                         ibool,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
!                         rho_k_m2,mu_k_m2,kappa_k_m2, & !lucas, CTD-SEM
!                         GPU_MODE,it,NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS, &
!                         !no_backward_displ_buffer_lucas_tmp,no_backward_displ_buffer_lucas, & ! lucas, this line needs to open when du* is computed
!                         ibool, &
!                         b_accel_elastic_m2 !lucas, CTD-SEM!

!  use specfem_par_gpu, only: Mesh_pointer,deltatf

!  use specfem_par, only: c11_k,c13_k,c15_k,c33_k,c35_k,c55_k,ispec_is_anisotropic, &
!                         rho_kl,c11_kl,c13_kl,c15_kl,c33_kl,c35_kl,c55_kl

!  implicit none

  !local variables
!  integer :: i,j,k,ispec,iglob
!  real(kind=CUSTOM_REAL) :: dux_dxi_m2,dux_dgamma_m2,duz_dxi_m2,duz_dgamma_m2 !lucas, CTD-SEM
!  real(kind=CUSTOM_REAL) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
!  real(kind=CUSTOM_REAL) :: dux_dxl_m2,dux_dzl_m2,duz_dxl_m2,duz_dzl_m2  !lucas, CTD-SEM
!  real(kind=CUSTOM_REAL) :: b_dux_dxi_m2,b_dux_dgamma_m2,b_duz_dxi_m2,b_duz_dgamma_m2 !lucas, CTD-SEM
!  real(kind=CUSTOM_REAL) :: b_dux_dxl, b_dux_dzl,b_duz_dxl,b_duz_dzl
!  real(kind=CUSTOM_REAL) :: b_dux_dxl_m2, b_dux_dzl_m2,b_duz_dxl_m2,b_duz_dzl_m2 !lucas, CTD-SEM
!  real(kind=CUSTOM_REAL) :: dsxx_m2,dsxz_m2,dszz_m2 !lucas, CTD-SEM
!  real(kind=CUSTOM_REAL) :: b_dsxx_m2,b_dsxz_m2,b_dszz_m2 !lucas, CTD-SEM
!  real(kind=CUSTOM_REAL) :: rhol_m2,mul_m2,kappal_m2 !lucas, CTD-SEM
!  real(kind=CUSTOM_REAL) :: rhol,mul,kappal 

  ! Jacobian matrix and determinant
!  double precision :: xixl,xizl,gammaxl,gammazl

  ! 1.lucas, to get kappa_k, mu_k, and rho_k at each iglob
  ! elastic kernels
 
!  if (.not. GPU_MODE) then
    ! updates kernels on CPU
!    do ispec = 1,nspec
!      if (ispec_is_elastic(ispec)) then
!        do j = 1,NGLLZ; do i = 1,NGLLX
!          ! derivative along x and along z

!          dux_dxi_m2 = 0._CUSTOM_REAL
!          duz_dxi_m2 = 0._CUSTOM_REAL
!          dux_dgamma_m2 = 0._CUSTOM_REAL
!          duz_dgamma_m2 = 0._CUSTOM_REAL
!          b_dux_dxi_m2 = 0._CUSTOM_REAL
!          b_duz_dxi_m2 = 0._CUSTOM_REAL
!          b_dux_dgamma_m2 = 0._CUSTOM_REAL
!          b_duz_dgamma_m2 = 0._CUSTOM_REAL


          ! 1.1 lucas, to get dux_dxi, duz_dxi, dux_dgamma, duz_dgamma and also the back conterparts
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX ! lucas: in Hessian Ha=du.u*, use 1.1.a; in Hessan Hb=u.du*, use 1.1.b; in frechet, use 1.1.a.
            
!            dux_dxi_m2 = dux_dxi_m2 + displ_elastic_m2(1,ibool(k,j,ispec))*hprime_xx(i,k)
!            duz_dxi_m2 = duz_dxi_m2 + displ_elastic_m2(2,ibool(k,j,ispec))*hprime_xx(i,k)
!            dux_dgamma_m2 = dux_dgamma_m2 + displ_elastic_m2(1,ibool(i,k,ispec))*hprime_zz(j,k)
!            duz_dgamma_m2 = duz_dgamma_m2 + displ_elastic_m2(2,ibool(i,k,ispec))*hprime_zz(j,k)
             
!            b_dux_dxi_m2 = b_dux_dxi_m2 + b_displ_elastic_m2(1,ibool(k,j,ispec))*hprime_xx(i,k)
!            b_duz_dxi_m2 = b_duz_dxi_m2 + b_displ_elastic_m2(2,ibool(k,j,ispec))*hprime_xx(i,k)
!            b_dux_dgamma_m2 = b_dux_dgamma_m2 + b_displ_elastic_m2(1,ibool(i,k,ispec))*hprime_zz(j,k)
!            b_duz_dgamma_m2 = b_duz_dgamma_m2 + b_displ_elastic_m2(2,ibool(i,k,ispec))*hprime_zz(j,k)
           

!          enddo
!          ! 1.2 lucas, to get xixl, xizl, gammaxl, gammazl  
!          xixl = xix(i,j,ispec)
!          xizl = xiz(i,j,ispec)
!          gammaxl = gammax(i,j,ispec)
!          gammazl = gammaz(i,j,ispec)
          ! 1.3 lucas, to get dux_dxl, dux_dzl, duz_dxl, duz_dzl, and also the back counterparts
          ! derivatives of displacement

!          dux_dxl_m2 = dux_dxi_m2*xixl + dux_dgamma_m2*gammaxl
!          dux_dzl_m2 = dux_dxi_m2*xizl + dux_dgamma_m2*gammazl
!          duz_dxl_m2 = duz_dxi_m2*xixl + duz_dgamma_m2*gammaxl
!          duz_dzl_m2 = duz_dxi_m2*xizl + duz_dgamma_m2*gammazl

!          b_dux_dxl_m2 = b_dux_dxi_m2*xixl + b_dux_dgamma_m2*gammaxl
!          b_dux_dzl_m2 = b_dux_dxi_m2*xizl + b_dux_dgamma_m2*gammazl

!          b_duz_dxl_m2 = b_duz_dxi_m2*xixl + b_duz_dgamma_m2*gammaxl
!          b_duz_dzl_m2 = b_duz_dxi_m2*xizl + b_duz_dgamma_m2*gammazl

!          iglob = ibool(i,j,ispec)
          ! 1.4 lucas dsxx, dsxz, dszz, the back counterparts, and kappa_k, mu_k
          ! isotropic kernel contributions
!          if (P_SV) then
            ! P-SV waves

            !lucas added for approximate hessian--------------------------------------------------------- 
            !kappa_k(iglob) = (b_dsxx + b_dszz) *  (b_dsxx + b_dszz)
            !mu_k(iglob) = b_dsxx * b_dsxx + b_dszz * b_dszz + &
            !              2._CUSTOM_REAL * b_dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob)
            ! --------------------------------------------------------------------------------------------


!            dsxx_m2 =  dux_dxl_m2
!            dsxz_m2 = HALF * (duz_dxl_m2 + dux_dzl_m2)
!            dszz_m2 =  duz_dzl_m2

!            b_dsxx_m2 =  b_dux_dxl_m2
!            b_dsxz_m2 = HALF * (b_duz_dxl_m2 + b_dux_dzl_m2)
!            b_dszz_m2 =  b_duz_dzl_m2

!            kappa_k_m2(iglob) = (dsxx_m2 + dszz_m2) *  (b_dsxx_m2 + b_dszz_m2)
!            mu_k_m2(iglob) = dsxx_m2 * b_dsxx_m2 + dszz_m2 * b_dszz_m2 + &
!                          2._CUSTOM_REAL * dsxz_m2 * b_dsxz_m2 - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k_m2(iglob)

!          else
            ! SH (membrane) waves
!            mu_k_m2(iglob) = dux_dxl_m2 * b_dux_dxl_m2 + dux_dzl_m2 * b_dux_dzl_m2

            ! lucas added for approximate hessian-----------------------
            ! mu_k(iglob) = b_dux_dxl * b_dux_dxl + b_dux_dzl * b_dux_dzl
            ! ----------------------------------------------------------
!          endif

          ! Voigt kernels, e.g., see Sieminski, 2007a,b
!          if (ispec_is_anisotropic(ispec)) then
!            c11_k(iglob) = dux_dxl*b_dux_dxl
!            c13_k(iglob) = dux_dxl*b_duz_dzl + duz_dzl*b_dux_dxl
!            c15_k(iglob) = 2*(dux_dxl*HALF*(b_dux_dzl+b_duz_dxl)+&
!                           HALF*(dux_dzl+duz_dxl)*b_dux_dxl)
!            c33_k(iglob) = duz_dzl*b_duz_dzl
!            c35_k(iglob) = 2*(duz_dzl*HALF*(b_dux_dzl+b_duz_dxl)+&
!                           HALF*(dux_dzl+duz_dxl)*b_duz_dzl)
!            c55_k(iglob) = 4*HALF*(dux_dzl+duz_dxl)*HALF*(b_dux_dzl+b_duz_dxl)
!          endif
!        enddo;enddo
!      endif
!    enddo

!    do iglob = 1,nglob
       !rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) + accel_elastic(2,iglob)*b_displ_elastic(2,iglob) ! lucas, where 1 means x, and 2 means z component.
       !lucas added for approximate Hessian------------------------------------------------------------------------------ 
       !rho_k(iglob) =  b_veloc_elastic(1,iglob)*b_veloc_elastic(1,iglob) + b_veloc_elastic(2,iglob)*b_veloc_elastic(2,iglob)
       ! ----------------------------------------------------------------------------------------------------------------

!       rho_k_m2(iglob) =  displ_elastic_m2(1,iglob)*b_accel_elastic_m2(1,iglob) + & 
!                          displ_elastic_m2(2,iglob)*b_accel_elastic_m2(2,iglob) ! lucas, correct version
   
!    enddo

!  else
    ! updates kernels on GPU
!    call compute_kernels_elastic_cuda(Mesh_pointer,deltatf)
!  endif

! !2 .lucas, to get kernels (rho_kl, mu_kl, kappa_kl) at each iglob
! do ispec = 1, nspec
!    if (ispec_is_elastic(ispec)) then
      ! isotropic kernels
!      do j = 1, NGLLZ
!        do i = 1, NGLLX
!          iglob = ibool(i,j,ispec)

          ! for parameterization (rho,mu,kappa): "primary" kernels
          ! density kernel  
!          rho_kl_m2(i,j,ispec) = rho_kl_m2(i,j,ispec) - rho_k_m2(iglob) 
          ! shear modulus kernel
!          mu_kl_m2(i,j,ispec) =  mu_kl_m2(i,j,ispec) -  mu_k_m2(iglob)
          ! bulk modulus kernel
!          kappa_kl_m2(i,j,ispec) = kappa_kl_m2(i,j,ispec) -  kappa_k_m2(iglob)
          

!        enddo
!      enddo
      ! Voigt kernels, e.g., see Sieminski, 2007a,b
!      if (ispec_is_anisotropic(ispec)) then
!        do j = 1, NGLLZ
!          do i = 1, NGLLX
!            iglob = ibool(i,j,ispec)
!            rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rho_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
!            c11_kl(i,j,ispec) = c11_kl(i,j,ispec) - c11_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
!            c13_kl(i,j,ispec) = c13_kl(i,j,ispec) - c13_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
!            c15_kl(i,j,ispec) = c15_kl(i,j,ispec) - c15_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
!            c33_kl(i,j,ispec) = c33_kl(i,j,ispec) - c33_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
!            c35_kl(i,j,ispec) = c35_kl(i,j,ispec) - c35_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
!            c55_kl(i,j,ispec) = c55_kl(i,j,ispec) - c55_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
!            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + c11_kl(i,j,ispec) + &
!                                 c13_kl(i,j,ispec) + c15_kl(i,j,ispec) + c33_kl(i,j,ispec) + &
!                                 c35_kl(i,j,ispec) + c55_kl(i,j,ispec)
!          enddo
!        enddo
!      endif


!    endif
!  enddo

  ! 3.lucas, multiply delta and parameeter at the last step since A0dt+A1dt+A2dt + ...+ Andt=(A0+A1+A2+...+An)dt, i.e., any constant at each time step can be done after the sum.
  !   lucas, we compute the kernels at a specfic it, which means tha this if only be execute when the condition is satisfied.
  ! only at the last time step we multiply by delta and parameter value, it is not necessary to do it at each iteration

!   if (NSTEP - it == mod(NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS)) then  !lucas:NSTEP_BETWEEN_COMPUTE_KERNELS=1 set in Par_file. mod(NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS) =0
  ! lucas test, we save the result at 7400 step, not until the last step. 
  !if (NSTEP-it == 2600) then ! lucas: 7400=10000-2600
!    do ispec = 1, nspec
!      if (ispec_is_elastic(ispec)) then
        ! isotropic kernels
!        do j = 1, NGLLZ
!          do i = 1, NGLLX

!            if (.not. assign_external_model) then
!              rhol = density(1,kmato(ispec))
!              mul = poroelastcoef(2,1,kmato(ispec))
              !if (AXISYM) then ! ABAB !!
              !Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
!                kappal = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS * mul
              !else
              !  kappal = poroelastcoef(3,1,kmato(ispec)) - mul
              !endif
!            else
!              rhol_m2 = rhoext_m2(i,j,ispec)
!              mul_m2 = rhoext_m2(i,j,ispec)*vsext_m2(i,j,ispec)*vsext_m2(i,j,ispec)
              !if (AXISYM) then ! ABAB !!
              ! Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
!              kappal_m2 = rhoext_m2(i,j,ispec)*vpext_m2(i,j,ispec)*vpext_m2(i,j,ispec) - FOUR_THIRDS * mul_m2
              !else
              !  kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - mul
              !endif
!            endif

            ! for parameterization (rho,mu,kappa): "primary" kernels
            ! density kernel
!            rho_kl_m2(i,j,ispec) = rhol_m2 * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl_m2(i,j,ispec) ! lucas see tromp et al. 2005
            ! shear modulus kernel
!            mu_kl_m2(i,j,ispec) =  TWO * mul_m2 * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl_m2(i,j,ispec)
            ! bulk modulus kernel
!            kappa_kl_m2(i,j,ispec) = kappal_m2 * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl_m2(i,j,ispec)

            ! for parameterization (rho,beta,alpha):
            ! rho prime kernel
!            rhop_kl_m2(i,j,ispec) = rho_kl_m2(i,j,ispec) + kappa_kl_m2(i,j,ispec) + mu_kl_m2(i,j,ispec)
            ! Vs kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
!            beta_kl_m2(i,j,ispec) = TWO * (mu_kl_m2(i,j,ispec) - FOUR_THIRDS * mul_m2/kappal_m2 * kappa_kl_m2(i,j,ispec))
            ! Vp kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
!            alpha_kl_m2(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul_m2/kappal_m2) * kappa_kl_m2(i,j,ispec)
            ! for bulk velocity c parameterization (rho,bulk_c,beta):
            !bulk_c_kl(i,j,ispec) =  TWO * kappa_kl(i,j,ispec) ! for binary format, see save_adjoint_kernels.f90
            !bulk_beta_kl(i,j,ispec) =  TWO * mu_kl(i,j,ispec) ! for binary format, see save_adjoint_kernels.f90


!          enddo
!        enddo
!      endif ! elastic
!    enddo !nspec loop
!  endif ! it == NSTEP

!  end subroutine compute_kernels_el_m2



!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_ac()

! acoustic kernel calculations
! see e.g. Tromp et al. (2005)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO,HALF,TWO

  use specfem_par, only: nspec,ispec_is_acoustic,ibool,kappal_ac_global,rhol_ac_global, &
                         poroelastcoef,density,kmato,assign_external_model,rhoext,vpext,deltat, &
                         hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         potential_acoustic,b_potential_acoustic,potential_dot_dot_acoustic, &
                         accel_ac,b_displ_ac,NSTEP_BETWEEN_COMPUTE_KERNELS, &
                         rho_ac_kl,kappa_ac_kl,rhop_ac_kl,alpha_ac_kl,GPU_MODE

  use specfem_par_gpu, only: Mesh_pointer,deltatf

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,b_tempx1l,b_tempx2l
  double precision :: xixl,xizl,gammaxl,gammazl

  if (.not. GPU_MODE) then
    ! kernels on CPU
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. assign_external_model) then
              kappal_ac_global(iglob) = poroelastcoef(3,1,kmato(ispec))
              rhol_ac_global(iglob) = density(1,kmato(ispec))
            else
              kappal_ac_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec)
              rhol_ac_global(iglob)   = rhoext(i,j,ispec)
            endif

            ! calcul the displacement by computing the gradient of potential / rho
            ! and calcul the acceleration by computing the gradient of potential_dot_dot / rho
            tempx1l = ZERO
            tempx2l = ZERO
            b_tempx1l = ZERO
            b_tempx2l = ZERO
          !  bb_tempx1l = ZERO
          !  bb_tempx2l = ZERO
            do k = 1,NGLLX
              ! derivative along x
              !tempx1l = tempx1l + potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
              tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k) !!! YANGL
              b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
              ! derivative along z
              !tempx2l = tempx2l + potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
              tempx2l = tempx2l + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k) !!! YANGL
              b_tempx2l = b_tempx2l + b_potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
            enddo

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)

            ! derivatives of potential
            accel_ac(1,iglob) = (tempx1l*xixl + tempx2l*gammaxl) / rhol_ac_global(iglob)
            accel_ac(2,iglob) = (tempx1l*xizl + tempx2l*gammazl) / rhol_ac_global(iglob)
            b_displ_ac(1,iglob) = (b_tempx1l*xixl + b_tempx2l*gammaxl) / rhol_ac_global(iglob)
            b_displ_ac(2,iglob) = (b_tempx1l*xizl + b_tempx2l*gammazl) / rhol_ac_global(iglob)
          enddo !i = 1, NGLLX
        enddo !j = 1, NGLLZ
      endif
    enddo

    do ispec = 1,nspec
      if (ispec_is_acoustic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            ! YANGL
            !!!! old expression (from elastic kernels)
            !!!rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol_ac_global(iglob)  * &
            !!!      dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
            !!!kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal_ac_global(iglob) * &
            !!!      potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
            !!!      b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob)&
            !!!      * deltat
            !!!! new expression (from PDE-constrained optimization, coupling terms changed as well)
            rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + rhol_ac_global(iglob) * &
                                   dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * &
                                   (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
                                   !warning : the variable is named accel_ac but it is displ_ac that is computed
            kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) + kappal_ac_global(iglob) * &
                                     potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
                                     b_potential_acoustic(iglob)/kappal_ac_global(iglob) * &
                                     (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
            ! YANGL
            rhop_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + kappa_ac_kl(i,j,ispec)
            alpha_ac_kl(i,j,ispec) = TWO *  kappa_ac_kl(i,j,ispec)
          enddo
        enddo
      endif
    enddo
  else
    ! on GPU
    call compute_kernels_acoustic_cuda(Mesh_pointer,deltatf)
  endif

  end subroutine compute_kernels_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_po()

! kernel calculations
! see e.g. Morency et al. (2009)

  use constants, only: CUSTOM_REAL,FOUR_THIRDS,NGLLX,NGLLZ,TWO,HALF

  use specfem_par, only: nglob,nspec,ispec_is_poroelastic,ibool,deltat, &
                         kmato,permeability, &
                         accels_poroelastic,accelw_poroelastic,velocw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic, &
                         epsilondev_s,b_epsilondev_s, &
                         epsilondev_w,b_epsilondev_w, &
                         rhot_k,rhof_k,sm_k,eta_k,B_k,C_k, &
                         rhot_kl,rhof_kl,sm_kl,eta_kl,B_kl,C_kl,M_kl,M_k, &
                         mufr_kl,mufr_k,rhob_kl,rhofb_kl, &
                         mufrb_kl,phi_kl,rhobb_kl,rhofbb_kl,phib_kl,cpI_kl,cpII_kl,cs_kl,ratio_kl, &
                         GPU_MODE,NSTEP_BETWEEN_COMPUTE_KERNELS
  implicit none

  !local variables
  integer :: i,j,ispec,iglob
  real(kind=CUSTOM_REAL) :: rholb,dd1
  real(kind=CUSTOM_REAL) :: ratio
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz,dszx_xz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz,b_dszx_xz
  real(kind=CUSTOM_REAL) :: dwxx,dwzz,b_dwxx,b_dwzz

  ! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: B_biot
  double precision :: perm_xx
  double precision :: afactor,bfactor,cfactor
  double precision :: gamma1,gamma2,gamma3,gamma4
  double precision :: cpIsquare,cpIIsquare,cssquare

  integer :: material

  ! safety check
  if (GPU_MODE) call stop_the_code('Error poroelastic kernels not implemented on GPUs yet')

  ! kernel contributions on global nodes
  do iglob = 1,nglob
    rhot_k(iglob) = accels_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                    accels_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob)

    rhof_k(iglob) = accelw_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                    accelw_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob) + &
                    accels_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                    accels_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)

    sm_k(iglob)  = accelw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                   accelw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)

    eta_k(iglob) = velocw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                   velocw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
  enddo

  ! kernels on local nodes
  do ispec = 1, nspec
    if (ispec_is_poroelastic(ispec)) then

      ! gets poroelastic material
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      B_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr

      ! permeability
      material = kmato(ispec)
      perm_xx = permeability(1,material)

      ! Approximated velocities (no viscous dissipation)
      afactor = rho_bar - phi/tort*rho_f
      bfactor = H_biot + phi*rho_bar/(tort*rho_f)*M_biot - TWO*phi/tort*C_biot
      cfactor = phi/(tort*rho_f)*(H_biot*M_biot - C_biot*C_biot)

      cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cssquare = mu_fr/afactor

      ! Approximated ratio r = amplitude "w" field/amplitude "s" field (no viscous dissipation)
      ! used later for wavespeed kernels calculation, which are presently implemented for inviscid case,
      ! contrary to primary and density-normalized kernels, which are consistent with viscous fluid case.
      gamma1 = H_biot - phi/tort*C_biot
      gamma2 = C_biot - phi/tort*M_biot
      gamma3 = phi/tort*( M_biot*(afactor/rho_f + phi/tort) - C_biot)
      gamma4 = phi/tort*( C_biot*(afactor/rho_f + phi/tort) - H_biot)

      ratio = HALF*(gamma1 - gamma3)/gamma4 + HALF*sqrt((gamma1-gamma3)**2/gamma4**2 + 4.d0 * gamma2/gamma4)

      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)

          rhot_kl(i,j,ispec) = rhot_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_bar * rhot_k(iglob)
          rhof_kl(i,j,ispec) = rhof_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_f * rhof_k(iglob)
          sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_f*tort/phi * sm_k(iglob)

          !at the moment works with constant permeability
          eta_kl(i,j,ispec) = eta_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * eta_f/perm_xx * eta_k(iglob)

          ! for B_k & mufr_k
          dsxx = epsilondev_s(1,i,j,ispec) ! dux_dxl
          dszz = epsilondev_s(2,i,j,ispec) ! duz_dzl
          dsxz = epsilondev_s(3,i,j,ispec) ! dux_dzl
          dszx_xz = epsilondev_s(4,i,j,ispec) ! 0.5_CUSTOM_REAL * (duz_dxl + dux_dzl)

          b_dsxx = b_epsilondev_s(1,i,j,ispec) ! b_dux_dxl
          b_dszz = b_epsilondev_s(2,i,j,ispec) ! b_duz_dzl
          b_dsxz = b_epsilondev_s(3,i,j,ispec) ! b_dux_dzl
          b_dszx_xz = b_epsilondev_s(4,i,j,ispec) ! 0.5_CUSTOM_REAL * (b_duz_dxl + b_dux_dzl)

          B_k(iglob) = (dsxx + dszz) *  (b_dsxx + b_dszz) * (H_biot - FOUR_THIRDS * mu_fr)

          mufr_k(iglob) = (dsxx * b_dsxx + dszz * b_dszz + &
                          2._CUSTOM_REAL * dszx_xz * b_dszx_xz - &
                          1._CUSTOM_REAL/3._CUSTOM_REAL * (dsxx + dszz) * (b_dsxx + b_dszz) ) * mu_fr

          ! from older compute_forces_poro_solid ...
          !  iglob = ibool(i,j,ispec)
          !  dsxx =  dux_dxl
          !  dsxz = HALF * (duz_dxl + dux_dzl)
          !  dszz =  duz_dzl
          !
          !  b_dsxx =  b_dux_dxl
          !  b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
          !  b_dszz =  b_duz_dzl
          !
          !  B_k(iglob) = (dux_dxl + duz_dzl) *  (b_dux_dxl + b_duz_dzl) * (H_biot - FOUR_THIRDS * mu_fr)
          !  mufr_k(iglob) = (dsxx * b_dsxx + dszz * b_dszz + &
          !                  2._CUSTOM_REAL * dsxz * b_dsxz - &
          !                  1._CUSTOM_REAL/3._CUSTOM_REAL * (dux_dxl + duz_dzl) * (b_dux_dxl + b_duz_dzl) ) * mu_fr

          ! for C_k & M_k
          dsxx = epsilondev_w(1,i,j,ispec) ! dux_dxl
          dszz = epsilondev_w(2,i,j,ispec) ! duz_dzl
          dwxx = epsilondev_w(3,i,j,ispec) ! dwx_dxl
          dwzz = epsilondev_w(4,i,j,ispec) ! dwz_dzl

          b_dsxx = b_epsilondev_w(1,i,j,ispec) ! b_dux_dxl
          b_dszz = b_epsilondev_w(2,i,j,ispec) ! b_duz_dzl
          b_dwxx = b_epsilondev_w(3,i,j,ispec) ! b_dwx_dxl
          b_dwzz = b_epsilondev_w(4,i,j,ispec) ! b_dwz_dzl

          C_k(iglob) =  ( (dsxx + dszz)*(b_dwxx + b_dwzz) + (dwxx + dwzz)*(b_dsxx + b_dszz) ) * C_biot
          M_k(iglob) = (dwxx + dwzz)*(b_dwxx + b_dwzz) * M_biot

          ! from older compute_forces_poro_fluid ...
          !C_k(iglob) =  ((dux_dxl + duz_dzl) *  (b_dwx_dxl + b_dwz_dzl) + &
          !                (dwx_dxl + dwz_dzl) *  (b_dux_dxl + b_duz_dzl)) * C_biot
          !M_k(iglob) = (dwx_dxl + dwz_dzl) *  (b_dwx_dxl + b_dwz_dzl) * M_biot

          B_kl(i,j,ispec) = B_kl(i,j,ispec) - (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * B_k(iglob)
          C_kl(i,j,ispec) = C_kl(i,j,ispec) - (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * C_k(iglob)
          M_kl(i,j,ispec) = M_kl(i,j,ispec) - (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * M_k(iglob)

          mufr_kl(i,j,ispec) = mufr_kl(i,j,ispec) - TWO * (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * mufr_k(iglob)

          ! density kernels
          rholb = rho_bar - phi*rho_f/tort
          rhob_kl(i,j,ispec) = rhot_kl(i,j,ispec) + B_kl(i,j,ispec) + mufr_kl(i,j,ispec)
          rhofb_kl(i,j,ispec) = rhof_kl(i,j,ispec) + C_kl(i,j,ispec) + M_kl(i,j,ispec) + sm_kl(i,j,ispec)

          mufrb_kl(i,j,ispec) = mufr_kl(i,j,ispec)
          phi_kl(i,j,ispec) = - sm_kl(i,j,ispec) - M_kl(i,j,ispec)

          ! wave speed kernels
          dd1 = (1._CUSTOM_REAL+rholb/rho_f)*ratio**2 + 2._CUSTOM_REAL*ratio + tort/phi

          rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) &
                - phi*rho_f/(tort*B_biot) * &
                  (cpIIsquare + (cpIsquare - cpIIsquare)*( (phi / &
                  tort*ratio +1._CUSTOM_REAL)/dd1 + &
                  (rho_bar**2*ratio**2/rho_f**2*(phi / tort*ratio+1._CUSTOM_REAL)*(phi/tort*ratio + &
                  phi/tort * &
                  (1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL) )/dd1**2 ) - &
                  FOUR_THIRDS*cssquare ) &
                  * B_kl(i,j,ispec) &
                - rho_bar*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                  (phi/tort*ratio + &
                  1._CUSTOM_REAL)**2/dd1**2*M_kl(i,j,ispec) + &
                  rho_bar*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                  (phi/tort*ratio+1._CUSTOM_REAL)/dd1 - &
                  phi*ratio/tort*(phi / tort*ratio+1._CUSTOM_REAL)*&
                  (1._CUSTOM_REAL+rho_bar*ratio/rho_f)/dd1**2) &
                  * C_kl(i,j,ispec) &
                + phi*rho_f*cssquare / (tort*mu_fr) &
                  * mufrb_kl(i,j,ispec)

          rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) &
                + phi*rho_f/(tort*B_biot) * (cpIIsquare + (cpIsquare - cpIIsquare)*( (phi/ &
                  tort*ratio +1._CUSTOM_REAL)/dd1+&
                  (rho_bar**2*ratio**2/rho_f**2*(phi/tort*ratio+1)*(phi/tort*ratio+ &
                  phi/tort*(1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL) )/dd1**2 )- &
                  FOUR_THIRDS*cssquare ) &
                  * B_kl(i,j,ispec) &
                + rho_bar*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                  (phi/tort*ratio + 1._CUSTOM_REAL)**2/dd1**2 &
                  * M_kl(i,j,ispec) &
                - rho_bar*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                  (phi/tort*ratio+1._CUSTOM_REAL)/dd1 - &
                  phi*ratio/tort*(phi/tort*ratio+1._CUSTOM_REAL)*&
                  (1._CUSTOM_REAL+rho_bar*ratio/rho_f)/dd1**2) &
                  * C_kl(i,j,ispec) &
                - phi*rho_f*cssquare/(tort*mu_fr) &
                  * mufrb_kl(i,j,ispec)

          phib_kl(i,j,ispec) = phi_kl(i,j,ispec) &
                - phi*rho_bar/(tort*B_biot) * ( cpIsquare - rho_f/rho_bar*cpIIsquare- &
                  (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phi/tort + (1._CUSTOM_REAL+rho_f/rho_bar)* &
                  (TWO*ratio*phi/tort+1._CUSTOM_REAL))/dd1 + (phi/tort*ratio+1._CUSTOM_REAL)*(phi*&
                  ratio/tort+phi/tort*(1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL)*((1._CUSTOM_REAL+ &
                  rho_bar/rho_f-TWO*phi/tort)*ratio**2+TWO*ratio)/dd1**2 ) - &
                  FOUR_THIRDS*rho_f*cssquare/rho_bar ) &
                  * B_kl(i,j,ispec) &
                + rho_f/M_biot * (cpIsquare-cpIIsquare) &
                  *( TWO*ratio*(phi/tort*ratio+1._CUSTOM_REAL)/dd1 - &
                    (phi/tort*ratio+1._CUSTOM_REAL)**2 &
                    *((1._CUSTOM_REAL+rho_bar/rho_f-TWO*phi/tort)*ratio**2+TWO*ratio)/dd1**2) &
                  * M_kl(i,j,ispec) &
                + phi*rho_f/(tort*C_biot)* (cpIsquare-cpIIsquare)*ratio* (&
                  (1._CUSTOM_REAL+rho_f/rho_bar*ratio)/dd1 - (phi/tort*ratio+1._CUSTOM_REAL)* &
                  (1._CUSTOM_REAL+rho_bar/rho_f*ratio)*((1._CUSTOM_REAL+rho_bar/rho_f-TWO*phi/tort)*ratio+TWO)/dd1**2 ) &
                  * C_kl(i,j,ispec) &
                - phi*rho_f*cssquare /(tort*mu_fr) &
                  * mufrb_kl(i,j,ispec)

          ! wavespeed kernels
          cpI_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rho_bar*( &
                  1._CUSTOM_REAL-phi/tort + (phi/tort*ratio+ 1._CUSTOM_REAL)*(phi/tort*&
                  ratio+phi/tort* (1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL)/dd1 ) &
                  * B_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIsquare*rho_f*tort/(phi*M_biot) *&
                  (phi/tort*ratio+1._CUSTOM_REAL)**2/dd1 &
                  * M_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIsquare*rho_f/C_biot * &
                  (phi/tort*ratio+1._CUSTOM_REAL)* (1._CUSTOM_REAL+rho_bar/rho_f*ratio)/dd1 &
                  * C_kl(i,j,ispec)
          cpII_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIIsquare*rho_bar/B_biot * (&
                  phi*rho_f/(tort*rho_bar) - (phi/tort*ratio+ 1._CUSTOM_REAL)*(phi/tort*ratio+phi/tort* &
                  (1._CUSTOM_REAL+rho_f/rho_bar)-&
                  1._CUSTOM_REAL)/dd1  ) &
                  * B_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIIsquare*rho_f*tort/(phi*M_biot) * (&
                  1._CUSTOM_REAL - (phi/tort*ratio+ 1._CUSTOM_REAL)**2/dd1  ) &
                  * M_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIIsquare*rho_f/C_biot * (&
                  1._CUSTOM_REAL - (phi/tort*ratio+ 1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                  rho_bar/rho_f*ratio)/dd1  ) &
                  * C_kl(i,j,ispec)

          cs_kl(i,j,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare* rho_bar/B_biot &
                  *(1._CUSTOM_REAL-phi*rho_f/(tort*rho_bar)) &
                  * B_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*(rho_bar-rho_f*phi/tort)/mu_fr*cssquare &
                  * mufrb_kl(i,j,ispec)

          ratio_kl(i,j,ispec) = ratio*rho_bar*phi/(tort*B_biot) * (cpIsquare-cpIIsquare) &
                  * (phi/tort*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rho_f/rho_bar)/dd1 - (phi/tort*ratio+1._CUSTOM_REAL)*&
                    (phi/tort*ratio+phi/tort*( 1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                      1._CUSTOM_REAL+rho_bar/rho_f-phi/tort) + 2._CUSTOM_REAL)/dd1**2  ) &
                  * B_kl(i,j,ispec) &
                + ratio*rho_f*tort/(phi*M_biot)*(cpIsquare-cpIIsquare) * 2._CUSTOM_REAL*phi/tort &
                  * ( (phi/tort*ratio+1._CUSTOM_REAL)/dd1 - (phi/tort*ratio+1._CUSTOM_REAL)**2 &
                      * ((1._CUSTOM_REAL+rho_bar/rho_f-phi/tort)*ratio + 1._CUSTOM_REAL)/dd1**2 ) &
                  * M_kl(i,j,ispec) &
                + ratio*rho_f/C_biot*(cpIsquare-cpIIsquare) &
                  * ( (2._CUSTOM_REAL*phi*rho_bar*ratio/(tort*rho_f)+phi/tort+rho_bar/rho_f)/dd1 - &
                       2._CUSTOM_REAL*phi/tort*(phi/tort*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+rho_bar/rho_f*ratio) &
                      *((1._CUSTOM_REAL + rho_bar/rho_f - phi/tort)*ratio+1._CUSTOM_REAL)/dd1**2 ) &
                  * C_kl(i,j,ispec)
        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_po

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_Hessian() !lucas, used for Yang Luo's P1 and P2

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nglob,nspec,ibool,ispec_is_acoustic,ispec_is_elastic, &
                         any_elastic,any_acoustic, &
                         accel_elastic,b_accel_elastic,accel_ac,b_accel_ac, &
                         rhorho_el_Hessian_final1,rhorho_el_Hessian_final2, &
                         rhorho_ac_Hessian_final1,rhorho_ac_Hessian_final2, &
                         deltat,GPU_MODE,NSTEP_BETWEEN_COMPUTE_KERNELS

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  !local variables
  real(kind=CUSTOM_REAL), dimension(nglob) :: rhorho_el_Hessian_temp1, rhorho_el_Hessian_temp2
  integer :: i,j,ispec,iglob


  if (.not. GPU_MODE) then
    ! elastic domains
    if (any_elastic) then
      ! approximate Hessians
      ! pre-computes contributions on global points
      do iglob = 1,nglob
        rhorho_el_Hessian_temp1(iglob) = b_accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &  ! P1=a times a
                                         b_accel_elastic(2,iglob)*b_accel_elastic(2,iglob)
        rhorho_el_Hessian_temp2(iglob) = accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &   !P2=a* times a
                                         accel_elastic(2,iglob)*b_accel_elastic(2,iglob)
      enddo

      ! on local GLL basis
      do ispec = 1, nspec
        if (ispec_is_elastic(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              rhorho_el_Hessian_final1(i,j,ispec) = rhorho_el_Hessian_final1(i,j,ispec) + &
                                                    rhorho_el_Hessian_temp1(iglob) * (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
              rhorho_el_Hessian_final2(i,j,ispec) = rhorho_el_Hessian_final2(i,j,ispec) + &
                                                    rhorho_el_Hessian_temp2(iglob) * (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
            enddo
          enddo
        endif
      enddo
    endif

    ! acoustic domains
    if (any_acoustic) then
      ! on local GLL basis
      do ispec = 1,nspec
        if (ispec_is_acoustic(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              rhorho_ac_Hessian_final1(i,j,ispec) = rhorho_ac_Hessian_final1(i,j,ispec) + &
                                                    dot_product(accel_ac(:,iglob),accel_ac(:,iglob)) * &
                                                    (deltat* NSTEP_BETWEEN_COMPUTE_KERNELS)
              rhorho_ac_Hessian_final2(i,j,ispec) = rhorho_ac_Hessian_final2(i,j,ispec) + &
                                                    dot_product(accel_ac(:,iglob),b_accel_ac(:,iglob)) * &
                                                    (deltat* NSTEP_BETWEEN_COMPUTE_KERNELS)
            enddo
          enddo
        endif
      enddo
    endif

  else
    ! on GPU
    ! computes contribution to density and bulk modulus kernel
    call compute_kernels_hess_cuda(Mesh_pointer,any_elastic,any_acoustic)
  endif

  end subroutine compute_kernels_Hessian

