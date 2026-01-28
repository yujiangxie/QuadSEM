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
!========================================================================

  subroutine compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old,dux_dxl_old,duz_dzl_old, &
                                         dux_dzl_plus_duz_dxl_old,PML_BOUNDARY_CONDITIONS,e1,e11,e13,iphase)

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    ONE,TWO,PI,TINYVAL,FOUR_THIRDS,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: nglob,assign_external_model,P_SV, &
                         ATTENUATION_VISCOELASTIC,nspec_ATT_el,N_SLS, &
                         ibool,kmato,ispec_is_elastic, &
                         poroelastcoef,xix,xiz,gammax,gammaz, &
                         jacobian,vpext,vsext,rhoext, &
                         QKappa_attenuation,Qmu_attenuation,QKappa_attenuationext,Qmu_attenuationext, &
                         c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,c22ext, &
                         ispec_is_anisotropic,anisotropy, &
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         it,coord,iglob_is_forced

  ! overlapping communication
  use specfem_par, only: nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  ! AXISYM
  use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic,veloc_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el,N_SLS),intent(inout) :: e1,e11,e13

  ! CPML coefficients and memory variables
  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  integer,intent(in) :: iphase

  !---
  !--- local variables
  !---
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: dummy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: tempx3
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: sigma_thetatheta
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempz1l,tempz2l
  real(kind=CUSTOM_REAL) :: theta,ct,st
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zy,sigma_zz,sigma_zx
  real(kind=CUSTOM_REAL) :: sigma_xx_prime,sigma_xz_prime,sigma_zz_prime,sigma_zx_prime
  real(kind=CUSTOM_REAL) :: xxi


  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum,e13_sum

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
    lambdaplus2mu_unrelaxed_elastic,cpl,csl,rhol,lambdalplusmul_unrelaxed_elastic

  ! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25,c22

  ! CPML coefficients and memory variables
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: accel_elastic_PML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  !additional variables for CPML implementation in viscoelastic simulation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl,mu_duz_dzl, &
                                                    mu_dux_dzl,mu_duz_dxl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
                                                    mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl,mu_pml_duz_dxl
  double precision :: qkappal,qmul

  integer :: num_elements,ispec_p

  ! this to avoid a warning at execution time about an undefined variable being used
  ! for the SH component in the case of a P-SV calculation, and vice versa
  ! P_SV-case
  sigma_xx = 0._CUSTOM_REAL
  sigma_xz = 0._CUSTOM_REAL
  sigma_zz = 0._CUSTOM_REAL
  sigma_zx = 0._CUSTOM_REAL
  ! SH-case
  sigma_xy = 0._CUSTOM_REAL
  sigma_zy = 0._CUSTOM_REAL

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  ! loop over spectral elements, lucas: to the end of the routine
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! only for elastic spectral elements
    if (.not. ispec_is_elastic(ispec)) cycle
    ! 1.lucas
    if (ATTENUATION_VISCOELASTIC) then
      if (.not. assign_external_model) then
        qkappal = QKappa_attenuation(kmato(ispec))
        qmul = Qmu_attenuation(kmato(ispec))
      else
        qkappal = maxval(QKappa_attenuationext(:,:,ispec))
        qmul =  maxval(Qmu_attenuationext(:,:,ispec))
      endif
    endif
    ! 2.lucas
    ! gets local displacement for element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        dummy_loc(1,i,j) = displ_elastic(1,iglob)
        dummy_loc(2,i,j) = displ_elastic(2,iglob)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
      enddo
    enddo
    !3.lucas
    ! first double loop over GLL points to compute and store gradients
    call mxm_4comp_singleA(dux_dxi,duz_dxi,dux_dgamma,duz_dgamma,dummy_loc,hprime_xx,hprime_zz) ! lucas: mxm_4comp_singleA(x1,x2,z1,z2,A,B,C), where out: x1,x2,z1,z2.
    !4.lucas
    ! AXISYM case overwrites dux_dxi and duz_dxi
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! derivative along x and along z
            dux_dxi(i,j) = 0._CUSTOM_REAL
            duz_dxi(i,j) = 0._CUSTOM_REAL
            do k = 1,NGLJ
              dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
              duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprimeBar_xx(i,k)
            enddo
          enddo
        enddo
      endif
    endif
    ! 5.lucas
    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)

        ! derivatives of displacement
        dux_dxl(i,j) = dux_dxi(i,j)*xixl + dux_dgamma(i,j)*gammaxl
        dux_dzl(i,j) = dux_dxi(i,j)*xizl + dux_dgamma(i,j)*gammazl

        duz_dxl(i,j) = duz_dxi(i,j)*xixl + duz_dgamma(i,j)*gammaxl
        duz_dzl(i,j) = duz_dxi(i,j)*xizl + duz_dgamma(i,j)*gammazl
      enddo
    enddo
    !6.lucas
    ! AXISYM case overwrite duz_dxl
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! d_uz/dr=0 on the axis
        ! i == 1
        do j = 1,NGLLZ
          duz_dxl(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif
    !7.lucas
    if (PML_BOUNDARY_CONDITIONS) then !viscoelastic PML
      if (ispec_is_PML(ispec)) then
        if (ATTENUATION_VISCOELASTIC) then
          if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then

            call pml_compute_memory_variables_viscoelastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                           kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
                                                           mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl, &
                                                           mu_pml_duz_dxl,kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl, &
                                                           mu_duz_dzl,mu_dux_dzl,mu_duz_dxl)
          else
            call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                     dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
                                                     PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
                                                     PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
          endif
        else
          call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                   dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
                                                   PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
                                                   PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
        endif
      endif ! ispec_is_PML
    endif ! PML_BOUNDARY_CONDITIONS
    !8.lucas
    ! AXISYM case overwrite duz_dxl and duz_dxl_prime
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! d_uz/dr=0 on the axis
        !i == 1
        do j = 1,NGLLZ
          duz_dxl(1,j) = 0._CUSTOM_REAL
          duz_dxl_prime(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif
    !9.lucas
    ! get unrelaxed elastic parameters of current spectral element
    lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
    lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))

    lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
    !10.lucas
    ! first double loop to compute gradient
    do j = 1,NGLLZ
      do i = 1,NGLLX
        !--- if external medium, get elastic parameters of current grid point
        ! 10.1 lucas
        if (assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          rhol = rhoext(i,j,ispec)

          mul_unrelaxed_elastic = rhol*csl*csl
          lambdal_unrelaxed_elastic = rhol*cpl*cpl - 2._CUSTOM_REAL * mul_unrelaxed_elastic
          lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic
          lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
        endif
        ! 10.2 lucas
        ! compute stress tensor (include attenuation or anisotropy if needed)
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            if (i == 1) then
              ! First GLJ point
              sigma_xx = 0._CUSTOM_REAL
              sigma_zz = 0._CUSTOM_REAL
              sigma_thetatheta(i,j) = 0._CUSTOM_REAL
              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              r_xiplus1(i,j) = xxi
              do k = 1,NGLJ
                sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
              enddo
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                        + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + sigma_xx/xxi)
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                       + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + sigma_zz/xxi)
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                 + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta(i,j)/xxi
            else
              ! Not first GLJ point
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                         + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                         + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                      + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
            endif
          else
            ! Not on the axis
            sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                       + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
            sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                       + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
            sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                    + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
          endif
        else ! Not axisym
          if (P_SV) then
            ! P_SV case
            sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * duz_dzl(i,j)
            sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * dux_dxl(i,j)
            sigma_xz = mul_unrelaxed_elastic * (duz_dxl(i,j) + dux_dzl(i,j))
            sigma_zx = sigma_xz
          else
            ! SH-case
            sigma_xy = mul_unrelaxed_elastic * dux_dxl(i,j)
            sigma_zy = mul_unrelaxed_elastic * dux_dzl(i,j)
          endif
        endif ! AXISYM

        ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
        ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below.

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

        ! 10.3 lucas
        if (ATTENUATION_VISCOELASTIC) then
          ! This routine updates the memory variables and with the current grad(displ)
          ! and gets the attenuation contribution (in e*_sum variables)
          call compute_attenuation_viscoelastic(e1,e11,e13,dux_dxl(i,j),dux_dzl(i,j),duz_dxl(i,j),duz_dzl(i,j), &
                                                dux_dxl_old,duz_dzl_old, &
                                                dux_dzl_plus_duz_dxl_old,PML_BOUNDARY_CONDITIONS,i,j,ispec, &
                                                e1_sum,e11_sum,e13_sum)

! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
          sigma_xx = sigma_xx + lambdalplusmul_unrelaxed_elastic * e1_sum + TWO * mul_unrelaxed_elastic * e11_sum
          sigma_xz = sigma_xz + mul_unrelaxed_elastic * e13_sum
          sigma_zz = sigma_zz + lambdalplusmul_unrelaxed_elastic * e1_sum - TWO * mul_unrelaxed_elastic * e11_sum
          sigma_zx = sigma_xz
        endif ! ATTENUATION_VISCOELASTIC
        ! 10.4 lucas
        if (PML_BOUNDARY_CONDITIONS) then
          if (ATTENUATION_VISCOELASTIC) then
            if (ispec_is_PML(ispec)) then
               if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then
                 sigma_xx = (lambdalplusmul_unrelaxed_elastic * kappa_pml_dux_dxl(i,j) + mul_unrelaxed_elastic * &
                            mu_pml_dux_dxl(i,j)) + (lambdalplusmul_unrelaxed_elastic * kappa_duz_dzl(i,j) - &
                            mul_unrelaxed_elastic * mu_duz_dzl(i,j))
                 sigma_xz = (mul_unrelaxed_elastic * mu_pml_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_dux_dzl(i,j))
                 sigma_zx = (mul_unrelaxed_elastic * mu_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_pml_dux_dzl(i,j))
                 sigma_zz = (lambdalplusmul_unrelaxed_elastic * kappa_dux_dxl(i,j) - mul_unrelaxed_elastic * &
                            mu_dux_dxl(i,j))+ (lambdalplusmul_unrelaxed_elastic * kappa_pml_duz_dzl(i,j) + &
                            mul_unrelaxed_elastic * mu_pml_duz_dzl(i,j))
               else
                 sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * PML_duz_dzl(i,j)
                 sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * PML_dux_dxl(i,j)
                 sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
                 sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
               endif
            endif
          else ! No ATTENUATION
            if (ispec_is_PML(ispec) .and. nspec_PML > 0) then
              if (ROTATE_PML_ACTIVATE) then
                theta = - ROTATE_PML_ANGLE/180._CUSTOM_REAL*PI
                if (it == 1) write(*,*) theta,ROTATE_PML_ACTIVATE,cos(theta),sin(theta)
                ct = cos(theta)
                st = sin(theta)
                sigma_xx_prime = lambdaplus2mu_unrelaxed_elastic * (ct**2*dux_dxl(i,j) + ct*st*duz_dxl(i,j) + &
                                                                    ct*st*dux_dzl(i,j) + st**2*duz_dzl(i,j)) &
                                 + lambdal_unrelaxed_elastic*(st**2*PML_dux_dxl(i,j) - ct*st*PML_duz_dxl(i,j) - &
                                                              ct*st*PML_dux_dzl(i,j) + ct**2*PML_duz_dzl(i,j))

                sigma_xz_prime = mul_unrelaxed_elastic * (-ct*st*dux_dxl(i,j) + ct**2*duz_dxl(i,j) - &
                                                          st**2*dux_dzl(i,j) + ct*st*duz_dzl(i,j)) &
                                 + mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) - st**2*PML_duz_dxl(i,j) + &
                                                            ct**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j))

                sigma_zx_prime = mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) + ct**2*PML_duz_dxl(i,j) - &
                                                          st**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j)) &
                                 + mul_unrelaxed_elastic * (-ct*st*dux_dxl_prime(i,j) - st**2*duz_dxl_prime(i,j) + &
                                                          ct**2*dux_dzl_prime(i,j) + ct*st*duz_dzl_prime(i,j))

                sigma_zz_prime = lambdaplus2mu_unrelaxed_elastic*(st**2*dux_dxl_prime(i,j) - ct*st*duz_dxl_prime(i,j) - &
                                                                  ct*st*dux_dzl_prime(i,j) + ct**2*duz_dzl_prime(i,j)) &
                                 + lambdal_unrelaxed_elastic*(ct**2*PML_dux_dxl(i,j) + ct*st*PML_duz_dxl(i,j) + &
                                                              ct*st*PML_dux_dzl(i,j) + st**2*PML_duz_dzl(i,j))

                sigma_xx = ct**2*sigma_xx_prime - ct*st*sigma_xz_prime - ct*st*sigma_zx_prime + st**2*sigma_zz_prime
                sigma_xz = ct*st*sigma_xx_prime + ct**2*sigma_xz_prime - st**2*sigma_zx_prime - ct*st*sigma_zz_prime
                sigma_zx = ct*st*sigma_xx_prime - st**2*sigma_xz_prime + ct**2*sigma_zx_prime - ct*st*sigma_zz_prime
                sigma_zz = st**2*sigma_xx_prime + ct*st*sigma_xz_prime + ct*st*sigma_zx_prime + ct**2*sigma_zz_prime
              else
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
                sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
                sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
              endif
            endif

          endif ! ATTENUATION
        endif ! PML_BOUNDARY_CONDITION
        ! 10.5 lucas
        ! full anisotropy
        if (ispec_is_anisotropic(ispec)) then
          if (assign_external_model) then
            c11 = c11ext(i,j,ispec)
            c13 = c13ext(i,j,ispec)
            c15 = c15ext(i,j,ispec)
            c33 = c33ext(i,j,ispec)
            c35 = c35ext(i,j,ispec)
            c55 = c55ext(i,j,ispec)
            c12 = c12ext(i,j,ispec)
            c23 = c23ext(i,j,ispec)
            c25 = c25ext(i,j,ispec)
            if (AXISYM) then
              c22 = c22ext(i,j,ispec) ! This variable is used for axisym simulations only
            endif
          else
            c11 = anisotropy(1,kmato(ispec))
            c13 = anisotropy(2,kmato(ispec))
            c15 = anisotropy(3,kmato(ispec))
            c33 = anisotropy(4,kmato(ispec))
            c35 = anisotropy(5,kmato(ispec))
            c55 = anisotropy(6,kmato(ispec))
            c12 = anisotropy(7,kmato(ispec))
            c23 = anisotropy(8,kmato(ispec))
            c25 = anisotropy(9,kmato(ispec))
            if (AXISYM) then
              c22 = anisotropy(10,kmato(ispec)) ! This variable is used for axisym simulations only
            endif
          endif

          ! implement anisotropy in 2D
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              if (i == 1) then
                ! first GLJ point, on the axis
                sigma_xx = 0._CUSTOM_REAL
                sigma_zz = 0._CUSTOM_REAL
                sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                r_xiplus1(i,j) = xxi
                do k = 1,NGLJ ! Compute the sum
                  sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                enddo
                sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*sigma_xx/xxi
                sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*sigma_zz/xxi
                sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*sigma_thetatheta(i,j)/xxi

              else
                ! not first GLJ point but not on the axis
                sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
                sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
                sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + &
                                        c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              endif
            else
              ! axisym but not on the axis
              sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
            endif
          else
            ! not AXISYM
            sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c15*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c35*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_zx = sigma_xz
          endif
        endif !end anisotropy

        ! 10.6 lucas
        ! weak formulation term based on stress tensor (non-symmetric form)
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        jacobianl = deriv(5,i,j)

        !! ABAB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
        ! tempx1(i,j) = w.J.F_{11}^{ij}
        ! tempz1(i,j) = w.J.F_{21}^{ij}
        ! tempx2(i,j) = w.J.F_{12}^{ij}
        ! tempz2(i,j) = w.J.F_{22}^{ij}
        ! 10.7 lucas
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            ! This is normal, we always add a contribution depending on the value on the axis
            ! i.e. we purposely sum something at point (i,j) with something at point (1,j)
            tempx3(i,j) = jacobian(1,j,ispec) * sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)

            ! not first GLJ point
            if (abs(coord(1,ibool(i,j,ispec))) > TINYVAL) then
              if (i == 1) then
                write(*,*) "Element number:",ispec
                call stop_the_code("Error: an axial element is rotated. The code should have been stopped before. Check that your &
                 &coordinates are greater than TINYVAL. Maybe you should also have a look to &
                 &doc/problematic_case_that_we_exclude_for_axisymmetric.pdf")
              endif
              tempx3(i,j) = tempx3(i,j) + wxglj(i) * jacobianl &
                            * sigma_thetatheta(i,j)/(xiglj(i)+ONE) ! this goes to accel_x
            endif

            tempx2(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
            tempx1(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
          else
            ! axisym but not on the axis
            tempx2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
            tempx1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
            tempx3(i,j) = wxgll(i) * jacobianl * sigma_thetatheta(i,j) ! this goes to accel_x
          endif
        else
          ! default (not axisym case)

          if (P_SV) then
            ! P_SV case
            tempx1(i,j) = jacobianl * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = jacobianl * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j) = jacobianl * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = jacobianl * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
          else
            ! SH-case
            tempx1(i,j) = jacobianl * (sigma_xy*xixl + sigma_zy*xizl) ! this goes to accel_x
            tempx2(i,j) = jacobianl * (sigma_xy*gammaxl + sigma_zy*gammazl) ! this goes to accel_x
            tempz1(i,j) = 0._CUSTOM_REAL
            tempz2(i,j) = 0._CUSTOM_REAL
          endif
        endif ! AXISYM
      enddo
    enddo  ! end of the loops on the collocation points i,j

    !11.lucas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update the displacement memory variable
    if (PML_BOUNDARY_CONDITIONS) then
      ! calculates contribution from each C-PML element to update acceleration
      call pml_compute_accel_contribution_elastic(ispec,nglob, &
                                                   dummy_loc,displ_elastic_old,veloc_elastic, &
                                                   accel_elastic_PML,r_xiplus1)
    endif
    !12.lucas
    !
    ! second double-loop over GLL to compute all the terms
    ! 12.1 lucas
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then

        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! along x direction and z direction
              ! and assemble the contributions
              ! we can merge the two loops because NGLLX == NGLLZ

              ! assembles the contributions
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              do k = 1,NGLLX
                tempx1l = tempx1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
                tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
                tempz1l = tempz1l + tempz1(k,j) * hprimeBarwglj_xx(k,i)
                tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxglj(i) * tempx2l)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxglj(i) * tempz2l)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
            endif
          enddo
        enddo

      else
        ! Axisym but not on the axis
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! assembles the contributions
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              do k = 1,NGLLX
                tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
                tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
                tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
                tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
            endif
          enddo
        enddo
      endif
    else ! 12.2 lucas
      ! default (not AXISYM case)
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (.not. iglob_is_forced(iglob)) then
            ! assembles the contributions
            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            do k = 1,NGLLX
              tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
              tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)
          endif
        enddo
      enddo
    endif ! AXISYM

    !13.lucas
    ! adds PML_BOUNDARY_CONDITIONS contribution
    if (PML_BOUNDARY_CONDITIONS) then
      if (ispec_is_PML(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - accel_elastic_PML(1,i,j)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) - accel_elastic_PML(2,i,j)
            endif
          enddo
        enddo
      endif
    endif

  enddo ! end of loop over all spectral elements

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_4comp_singleA(x1,x2,z1,z2,A,B,C)

! matrix x matrix multiplication, merging 4 loops for x1,x2 = A^t B^t and z1,z2 = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          duz_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!          duz_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprime_xx(i,k)
!            duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + dummy_loc(1,i,k)*hprime_zz(j,k)
!            duz_dgamma(i,j) = duz_dgamma(i,j) + dummy_loc(2,i,k)*hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NDIM,NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x1,x2,z1,z2

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: A
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x1(i,j) = A(1,1,j) * B(i,1) + A(1,2,j) * B(i,2) + A(1,3,j) * B(i,3) + A(1,4,j) * B(i,4) + A(1,5,j) * B(i,5)
        x2(i,j) = A(2,1,j) * B(i,1) + A(2,2,j) * B(i,2) + A(2,3,j) * B(i,3) + A(2,4,j) * B(i,4) + A(2,5,j) * B(i,5)

        z1(i,j) = A(1,i,1) * C(j,1) + A(1,i,2) * C(j,2) + A(1,i,3) * C(j,3) + A(1,i,4) * C(j,4) + A(1,i,5) * C(j,5)
        z2(i,j) = A(2,i,1) * C(j,1) + A(2,i,2) * C(j,2) + A(2,i,3) * C(j,3) + A(2,i,4) * C(j,4) + A(2,i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x1(i,j) = 0._CUSTOM_REAL
        x2(i,j) = 0._CUSTOM_REAL
        z1(i,j) = 0._CUSTOM_REAL
        z2(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x1(i,j) = x1(i,j) + A(1,k,j) * B(i,k)
          x2(i,j) = x2(i,j) + A(2,k,j) * B(i,k)

          z1(i,j) = z1(i,j) + A(1,i,k) * C(j,k)
          z2(i,j) = z2(i,j) + A(2,i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_4comp_singleA

  end subroutine compute_forces_viscoelastic

 !-------------lucas add it for backward, only for test, not used

  subroutine compute_forces_viscoelastic_backward(b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old,& 
                                                  b_dux_dxl_old,b_duz_dzl_old,b_dux_dzl_plus_duz_dxl_old,& 
                                                  PML_BOUNDARY_CONDITIONS,b_e1,b_e11,b_e13,iphase)

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    ONE,TWO,PI,TINYVAL,FOUR_THIRDS,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: nglob,assign_external_model,P_SV, &
                         ATTENUATION_VISCOELASTIC,nspec_ATT_el,N_SLS, &
                         ibool,kmato,ispec_is_elastic, &
                         poroelastcoef,xix,xiz,gammax,gammaz, &
                         jacobian,vpext,vsext,rhoext, &
                         QKappa_attenuation,Qmu_attenuation,QKappa_attenuationext,Qmu_attenuationext, &
                         c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,c22ext, &
                         ispec_is_anisotropic,anisotropy, &
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         it,coord,iglob_is_forced

  ! overlapping communication
  use specfem_par, only: nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  ! AXISYM
  use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: b_displ_elastic,b_veloc_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: b_accel_elastic

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: b_displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: b_dux_dxl_old,b_duz_dzl_old,&
                                                                               b_dux_dzl_plus_duz_dxl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el,N_SLS),intent(inout) :: b_e1,b_e11,b_e13

  ! CPML coefficients and memory variables
  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  integer,intent(in) :: iphase

  !---
  !--- local variables
  !---
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: dummy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: tempx3
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: sigma_thetatheta
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempz1l,tempz2l
  real(kind=CUSTOM_REAL) :: theta,ct,st
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zy,sigma_zz,sigma_zx
  real(kind=CUSTOM_REAL) :: sigma_xx_prime,sigma_xz_prime,sigma_zz_prime,sigma_zx_prime
  real(kind=CUSTOM_REAL) :: xxi


  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: b_e1_sum,b_e11_sum,b_e13_sum

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
    lambdaplus2mu_unrelaxed_elastic,cpl,csl,rhol,lambdalplusmul_unrelaxed_elastic

  ! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25,c22

  ! CPML coefficients and memory variables
  !real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: accel_elastic_PML, lucas
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  !additional variables for CPML implementation in viscoelastic simulation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl,mu_duz_dzl, &
                                                    mu_dux_dzl,mu_duz_dxl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
                                                    mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl,mu_pml_duz_dxl
  double precision :: qkappal,qmul

  integer :: num_elements,ispec_p

  ! this to avoid a warning at execution time about an undefined variable being used
  ! for the SH component in the case of a P-SV calculation, and vice versa
  ! P_SV-case
  sigma_xx = 0._CUSTOM_REAL
  sigma_xz = 0._CUSTOM_REAL
  sigma_zz = 0._CUSTOM_REAL
  sigma_zx = 0._CUSTOM_REAL
  ! SH-case
  sigma_xy = 0._CUSTOM_REAL
  sigma_zy = 0._CUSTOM_REAL

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  ! loop over spectral elements, lucas: to the end of the routine
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! only for elastic spectral elements
    if (.not. ispec_is_elastic(ispec)) cycle
    ! 1.lucas
    if (ATTENUATION_VISCOELASTIC) then
      if (.not. assign_external_model) then
        qkappal = QKappa_attenuation(kmato(ispec))
        qmul = Qmu_attenuation(kmato(ispec))
      else
        qkappal = maxval(QKappa_attenuationext(:,:,ispec))
        qmul =  maxval(Qmu_attenuationext(:,:,ispec))
      endif
    endif
    ! 2.lucas
    ! gets local displacement for element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        dummy_loc(1,i,j) = b_displ_elastic(1,iglob)
        dummy_loc(2,i,j) = b_displ_elastic(2,iglob)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
      enddo
    enddo
    !3.lucas
    ! first double loop over GLL points to compute and store gradients
    call mxm_4comp_singleA(dux_dxi,duz_dxi,dux_dgamma,duz_dgamma,dummy_loc,hprime_xx,hprime_zz) ! lucas: mxm_4comp_singleA(x1,x2,z1,z2,A,B,C), where out: x1,x2,z1,z2.
    !4.lucas
    ! AXISYM case overwrites dux_dxi and duz_dxi
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! derivative along x and along z
            dux_dxi(i,j) = 0._CUSTOM_REAL
            duz_dxi(i,j) = 0._CUSTOM_REAL
            do k = 1,NGLJ
              dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
              duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprimeBar_xx(i,k)
            enddo
          enddo
        enddo
      endif
    endif
    ! 5.lucas
    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)

        ! derivatives of displacement
        dux_dxl(i,j) = dux_dxi(i,j)*xixl + dux_dgamma(i,j)*gammaxl
        dux_dzl(i,j) = dux_dxi(i,j)*xizl + dux_dgamma(i,j)*gammazl

        duz_dxl(i,j) = duz_dxi(i,j)*xixl + duz_dgamma(i,j)*gammaxl
        duz_dzl(i,j) = duz_dxi(i,j)*xizl + duz_dgamma(i,j)*gammazl
      enddo
    enddo
    !6.lucas
    ! AXISYM case overwrite duz_dxl
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! d_uz/dr=0 on the axis
        ! i == 1
        do j = 1,NGLLZ
          duz_dxl(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif
    !7.lucas
    if (PML_BOUNDARY_CONDITIONS) then !viscoelastic PML
      if (ispec_is_PML(ispec)) then
        if (ATTENUATION_VISCOELASTIC) then
          if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then

            call pml_compute_memory_variables_viscoelastic(ispec,nglob,b_displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                           kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
                                                           mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl, &
                                                           mu_pml_duz_dxl,kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl, &
                                                           mu_duz_dzl,mu_dux_dzl,mu_duz_dxl)
          else
            call pml_compute_memory_variables_elastic(ispec,nglob,b_displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                     dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
                                                     PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
                                                     PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
          endif
        else
          call pml_compute_memory_variables_elastic(ispec,nglob,b_displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                   dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
                                                   PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
                                                   PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
        endif
      endif ! ispec_is_PML
    endif ! PML_BOUNDARY_CONDITIONS
    !8.lucas
    ! AXISYM case overwrite duz_dxl and duz_dxl_prime
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! d_uz/dr=0 on the axis
        !i == 1
        do j = 1,NGLLZ
          duz_dxl(1,j) = 0._CUSTOM_REAL
          duz_dxl_prime(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif
    !9.lucas
    ! get unrelaxed elastic parameters of current spectral element
    lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
    lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))

    lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
    !10.lucas
    ! first double loop to compute gradient
    do j = 1,NGLLZ
      do i = 1,NGLLX
        !--- if external medium, get elastic parameters of current grid point
        ! 10.1 lucas
        if (assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          rhol = rhoext(i,j,ispec)

          mul_unrelaxed_elastic = rhol*csl*csl
          lambdal_unrelaxed_elastic = rhol*cpl*cpl - 2._CUSTOM_REAL * mul_unrelaxed_elastic
          lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic
          lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
        endif
        ! 10.2 lucas
        ! compute stress tensor (include attenuation or anisotropy if needed)
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            if (i == 1) then
              ! First GLJ point
              sigma_xx = 0._CUSTOM_REAL
              sigma_zz = 0._CUSTOM_REAL
              sigma_thetatheta(i,j) = 0._CUSTOM_REAL
              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              r_xiplus1(i,j) = xxi
              do k = 1,NGLJ
                sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
              enddo
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                        + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + sigma_xx/xxi)
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                       + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + sigma_zz/xxi)
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                 + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta(i,j)/xxi
            else
              ! Not first GLJ point
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                         + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                         + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                      + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
            endif
          else
            ! Not on the axis
            sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                       + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
            sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                       + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
            sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                    + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
          endif
        else ! Not axisym
          if (P_SV) then
            ! P_SV case
            sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * duz_dzl(i,j)
            sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * dux_dxl(i,j)
            sigma_xz = mul_unrelaxed_elastic * (duz_dxl(i,j) + dux_dzl(i,j))
            sigma_zx = sigma_xz
          else
            ! SH-case
            sigma_xy = mul_unrelaxed_elastic * dux_dxl(i,j)
            sigma_zy = mul_unrelaxed_elastic * dux_dzl(i,j)
          endif
        endif ! AXISYM

        ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
        ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below.

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

        ! 10.3 lucas
        if (ATTENUATION_VISCOELASTIC) then
          ! This routine updates the memory variables and with the current grad(displ)
          ! and gets the attenuation contribution (in e*_sum variables)
          call compute_attenuation_viscoelastic(b_e1,b_e11,b_e13,dux_dxl(i,j),dux_dzl(i,j),duz_dxl(i,j),duz_dzl(i,j), &
                                                b_dux_dxl_old,b_duz_dzl_old, &
                                                b_dux_dzl_plus_duz_dxl_old,PML_BOUNDARY_CONDITIONS,i,j,ispec, &
                                                b_e1_sum,b_e11_sum,b_e13_sum)

! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
          sigma_xx = sigma_xx + lambdalplusmul_unrelaxed_elastic * b_e1_sum + TWO * mul_unrelaxed_elastic * b_e11_sum
          sigma_xz = sigma_xz + mul_unrelaxed_elastic * b_e13_sum
          sigma_zz = sigma_zz + lambdalplusmul_unrelaxed_elastic * b_e1_sum - TWO * mul_unrelaxed_elastic * b_e11_sum
          sigma_zx = sigma_xz
        endif ! ATTENUATION_VISCOELASTIC
        ! 10.4 lucas
        if (PML_BOUNDARY_CONDITIONS) then ! lucas, no considered yet
          if (ATTENUATION_VISCOELASTIC) then
            if (ispec_is_PML(ispec)) then
               if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then
                 sigma_xx = (lambdalplusmul_unrelaxed_elastic * kappa_pml_dux_dxl(i,j) + mul_unrelaxed_elastic * &
                            mu_pml_dux_dxl(i,j)) + (lambdalplusmul_unrelaxed_elastic * kappa_duz_dzl(i,j) - &
                            mul_unrelaxed_elastic * mu_duz_dzl(i,j))
                 sigma_xz = (mul_unrelaxed_elastic * mu_pml_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_dux_dzl(i,j))
                 sigma_zx = (mul_unrelaxed_elastic * mu_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_pml_dux_dzl(i,j))
                 sigma_zz = (lambdalplusmul_unrelaxed_elastic * kappa_dux_dxl(i,j) - mul_unrelaxed_elastic * &
                            mu_dux_dxl(i,j))+ (lambdalplusmul_unrelaxed_elastic * kappa_pml_duz_dzl(i,j) + &
                            mul_unrelaxed_elastic * mu_pml_duz_dzl(i,j))
               else
                 sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * PML_duz_dzl(i,j)
                 sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * PML_dux_dxl(i,j)
                 sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
                 sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
               endif
            endif
          else ! No ATTENUATION
            if (ispec_is_PML(ispec) .and. nspec_PML > 0) then
              if (ROTATE_PML_ACTIVATE) then
                theta = - ROTATE_PML_ANGLE/180._CUSTOM_REAL*PI
                if (it == 1) write(*,*) theta,ROTATE_PML_ACTIVATE,cos(theta),sin(theta)
                ct = cos(theta)
                st = sin(theta)
                sigma_xx_prime = lambdaplus2mu_unrelaxed_elastic * (ct**2*dux_dxl(i,j) + ct*st*duz_dxl(i,j) + &
                                                                    ct*st*dux_dzl(i,j) + st**2*duz_dzl(i,j)) &
                                 + lambdal_unrelaxed_elastic*(st**2*PML_dux_dxl(i,j) - ct*st*PML_duz_dxl(i,j) - &
                                                              ct*st*PML_dux_dzl(i,j) + ct**2*PML_duz_dzl(i,j))

                sigma_xz_prime = mul_unrelaxed_elastic * (-ct*st*dux_dxl(i,j) + ct**2*duz_dxl(i,j) - &
                                                          st**2*dux_dzl(i,j) + ct*st*duz_dzl(i,j)) &
                                 + mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) - st**2*PML_duz_dxl(i,j) + &
                                                            ct**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j))

                sigma_zx_prime = mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) + ct**2*PML_duz_dxl(i,j) - &
                                                          st**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j)) &
                                 + mul_unrelaxed_elastic * (-ct*st*dux_dxl_prime(i,j) - st**2*duz_dxl_prime(i,j) + &
                                                          ct**2*dux_dzl_prime(i,j) + ct*st*duz_dzl_prime(i,j))

                sigma_zz_prime = lambdaplus2mu_unrelaxed_elastic*(st**2*dux_dxl_prime(i,j) - ct*st*duz_dxl_prime(i,j) - &
                                                                  ct*st*dux_dzl_prime(i,j) + ct**2*duz_dzl_prime(i,j)) &
                                 + lambdal_unrelaxed_elastic*(ct**2*PML_dux_dxl(i,j) + ct*st*PML_duz_dxl(i,j) + &
                                                              ct*st*PML_dux_dzl(i,j) + st**2*PML_duz_dzl(i,j))

                sigma_xx = ct**2*sigma_xx_prime - ct*st*sigma_xz_prime - ct*st*sigma_zx_prime + st**2*sigma_zz_prime
                sigma_xz = ct*st*sigma_xx_prime + ct**2*sigma_xz_prime - st**2*sigma_zx_prime - ct*st*sigma_zz_prime
                sigma_zx = ct*st*sigma_xx_prime - st**2*sigma_xz_prime + ct**2*sigma_zx_prime - ct*st*sigma_zz_prime
                sigma_zz = st**2*sigma_xx_prime + ct*st*sigma_xz_prime + ct*st*sigma_zx_prime + ct**2*sigma_zz_prime
              else
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
                sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
                sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
              endif
            endif

          endif ! ATTENUATION
        endif ! PML_BOUNDARY_CONDITION
        ! 10.5 lucas
        ! full anisotropy
        if (ispec_is_anisotropic(ispec)) then
          if (assign_external_model) then
            c11 = c11ext(i,j,ispec)
            c13 = c13ext(i,j,ispec)
            c15 = c15ext(i,j,ispec)
            c33 = c33ext(i,j,ispec)
            c35 = c35ext(i,j,ispec)
            c55 = c55ext(i,j,ispec)
            c12 = c12ext(i,j,ispec)
            c23 = c23ext(i,j,ispec)
            c25 = c25ext(i,j,ispec)
            if (AXISYM) then
              c22 = c22ext(i,j,ispec) ! This variable is used for axisym simulations only
            endif
          else
            c11 = anisotropy(1,kmato(ispec))
            c13 = anisotropy(2,kmato(ispec))
            c15 = anisotropy(3,kmato(ispec))
            c33 = anisotropy(4,kmato(ispec))
            c35 = anisotropy(5,kmato(ispec))
            c55 = anisotropy(6,kmato(ispec))
            c12 = anisotropy(7,kmato(ispec))
            c23 = anisotropy(8,kmato(ispec))
            c25 = anisotropy(9,kmato(ispec))
            if (AXISYM) then
              c22 = anisotropy(10,kmato(ispec)) ! This variable is used for axisym simulations only
            endif
          endif

          ! implement anisotropy in 2D
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              if (i == 1) then
                ! first GLJ point, on the axis
                sigma_xx = 0._CUSTOM_REAL
                sigma_zz = 0._CUSTOM_REAL
                sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                r_xiplus1(i,j) = xxi
                do k = 1,NGLJ ! Compute the sum
                  sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                enddo
                sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*sigma_xx/xxi
                sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*sigma_zz/xxi
                sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*sigma_thetatheta(i,j)/xxi

              else
                ! not first GLJ point but not on the axis
                sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
                sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
                sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + &
                                        c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              endif
            else
              ! axisym but not on the axis
              sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
            endif
          else
            ! not AXISYM
            sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c15*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c35*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
            sigma_zx = sigma_xz
          endif
        endif
        ! 10.6 lucas
        ! weak formulation term based on stress tensor (non-symmetric form)
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        jacobianl = deriv(5,i,j)

        !! ABAB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
        ! tempx1(i,j) = w.J.F_{11}^{ij}
        ! tempz1(i,j) = w.J.F_{21}^{ij}
        ! tempx2(i,j) = w.J.F_{12}^{ij}
        ! tempz2(i,j) = w.J.F_{22}^{ij}
        ! 10.7 lucas
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            ! This is normal, we always add a contribution depending on the value on the axis
            ! i.e. we purposely sum something at point (i,j) with something at point (1,j)
            tempx3(i,j) = jacobian(1,j,ispec) * sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)

            ! not first GLJ point
            if (abs(coord(1,ibool(i,j,ispec))) > TINYVAL) then
              if (i == 1) then
                write(*,*) "Element number:",ispec
                call stop_the_code("Error: an axial element is rotated. The code should have been stopped before. Check that your &
                 &coordinates are greater than TINYVAL. Maybe you should also have a look to &
                 &doc/problematic_case_that_we_exclude_for_axisymmetric.pdf")
              endif
              tempx3(i,j) = tempx3(i,j) + wxglj(i) * jacobianl &
                            * sigma_thetatheta(i,j)/(xiglj(i)+ONE) ! this goes to accel_x
            endif

            tempx2(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
            tempx1(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
          else
            ! axisym but not on the axis
            tempx2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
            tempx1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
            tempx3(i,j) = wxgll(i) * jacobianl * sigma_thetatheta(i,j) ! this goes to accel_x
          endif
        else
          ! default (not axisym case)

          if (P_SV) then
            ! P_SV case
            tempx1(i,j) = jacobianl * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = jacobianl * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j) = jacobianl * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = jacobianl * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
          else
            ! SH-case
            tempx1(i,j) = jacobianl * (sigma_xy*xixl + sigma_zy*xizl) ! this goes to accel_x
            tempx2(i,j) = jacobianl * (sigma_xy*gammaxl + sigma_zy*gammazl) ! this goes to accel_x
            tempz1(i,j) = 0._CUSTOM_REAL
            tempz2(i,j) = 0._CUSTOM_REAL
          endif
        endif ! AXISYM
      enddo
    enddo  ! end of the loops on the collocation points i,j
    !11.lucas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update the displacement memory variable
   ! if (PML_BOUNDARY_CONDITIONS) then
   !   ! calculates contribution from each C-PML element to update acceleration
   !   call pml_compute_accel_contribution_elastic(ispec,nglob, &
   !                                                dummy_loc,displ_elastic_old,veloc_elastic, &
   !                                                accel_elastic_PML,r_xiplus1)
   ! endif
    !12.lucas
    !
    ! second double-loop over GLL to compute all the terms
    ! 12.1 lucas
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then

        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! along x direction and z direction
              ! and assemble the contributions
              ! we can merge the two loops because NGLLX == NGLLZ

              ! assembles the contributions
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              do k = 1,NGLLX
                tempx1l = tempx1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
                tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
                tempz1l = tempz1l + tempz1(k,j) * hprimeBarwglj_xx(k,i)
                tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxglj(i) * tempx2l)
              b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxglj(i) * tempz2l)

              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
            endif
          enddo
        enddo

      else
        ! Axisym but not on the axis
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! assembles the contributions
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              do k = 1,NGLLX
                tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
                tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
                tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
                tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
              b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)

              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
            endif
          enddo
        enddo
      endif
    else ! 12.2 lucas
      ! default (not AXISYM case)
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (.not. iglob_is_forced(iglob)) then
            ! assembles the contributions
            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            do k = 1,NGLLX
              tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
              tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
            b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)
          endif
        enddo
      enddo
    endif ! AXISYM

    !13.lucas
    ! adds PML_BOUNDARY_CONDITIONS contribution
  !  if (PML_BOUNDARY_CONDITIONS) then
  !    if (ispec_is_PML(ispec)) then
  !      do j = 1,NGLLZ
  !        do i = 1,NGLLX
  !          iglob = ibool(i,j,ispec)
  !          if (.not. iglob_is_forced(iglob)) then
  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - accel_elastic_PML(1,i,j)
  !            accel_elastic(2,iglob) = accel_elastic(2,iglob) - accel_elastic_PML(2,i,j)
  !          endif
  !        enddo
  !      enddo
  !    endif
  !  endif

  enddo ! end of loop over all spectral elements

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_4comp_singleA(x1,x2,z1,z2,A,B,C)

! matrix x matrix multiplication, merging 4 loops for x1,x2 = A^t B^t and z1,z2 = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          duz_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!          duz_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprime_xx(i,k)
!            duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + dummy_loc(1,i,k)*hprime_zz(j,k)
!            duz_dgamma(i,j) = duz_dgamma(i,j) + dummy_loc(2,i,k)*hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NDIM,NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x1,x2,z1,z2

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: A
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x1(i,j) = A(1,1,j) * B(i,1) + A(1,2,j) * B(i,2) + A(1,3,j) * B(i,3) + A(1,4,j) * B(i,4) + A(1,5,j) * B(i,5)
        x2(i,j) = A(2,1,j) * B(i,1) + A(2,2,j) * B(i,2) + A(2,3,j) * B(i,3) + A(2,4,j) * B(i,4) + A(2,5,j) * B(i,5)

        z1(i,j) = A(1,i,1) * C(j,1) + A(1,i,2) * C(j,2) + A(1,i,3) * C(j,3) + A(1,i,4) * C(j,4) + A(1,i,5) * C(j,5)
        z2(i,j) = A(2,i,1) * C(j,1) + A(2,i,2) * C(j,2) + A(2,i,3) * C(j,3) + A(2,i,4) * C(j,4) + A(2,i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x1(i,j) = 0._CUSTOM_REAL
        x2(i,j) = 0._CUSTOM_REAL
        z1(i,j) = 0._CUSTOM_REAL
        z2(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x1(i,j) = x1(i,j) + A(1,k,j) * B(i,k)
          x2(i,j) = x2(i,j) + A(2,k,j) * B(i,k)

          z1(i,j) = z1(i,j) + A(1,i,k) * C(j,k)
          z2(i,j) = z2(i,j) + A(2,i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_4comp_singleA

  end subroutine compute_forces_viscoelastic_backward !lucas, add it





!! ----------------------------------------------------------lucas, CTD-SEM method -------------------------------------------------------------
!! Yujiang Lucas Xie

  subroutine compute_forces_viscoelastic_m2(accel_elastic_m2,veloc_elastic_m2,displ_elastic_m2,displ_elastic_old_m2, & 
                                         dux_dxl_old_m2,duz_dzl_old_m2,dux_dzl_plus_duz_dxl_old_m2, & 
                                         PML_BOUNDARY_CONDITIONS,e1_m2,e11_m2,e13_m2,iphase) ! lucas, CTD-SEM
                                                                                     

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    ONE,TWO,PI,TINYVAL,FOUR_THIRDS,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: nglob,assign_external_model,P_SV, &
                         ATTENUATION_VISCOELASTIC,nspec_ATT_el,N_SLS, &
                         ibool,kmato,ispec_is_elastic, &
                         xix,xiz,gammax,gammaz, &
                         jacobian,vpext_m2,vsext_m2,rhoext_m2, & ! lucas, CTD-SEM
                         QKappa_attenuation,Qmu_attenuation,QKappa_attenuationext_m2,Qmu_attenuationext_m2, & ! lucas, CTD-SEM
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         iglob_is_forced

  ! overlapping communication
  use specfem_par, only: nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  ! AXISYM
 ! use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  ! PML arrays
 ! use specfem_par, only: nspec_PML,ispec_is_PML,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_m2,veloc_elastic_m2 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic_m2 !lucas, CTD-SEM

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_old_m2 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: dux_dxl_old_m2,duz_dzl_old_m2 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: dux_dzl_plus_duz_dxl_old_m2 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el,N_SLS),intent(inout) :: e1_m2,e11_m2,e13_m2 !lucas, CTD-SEM

  ! CPML coefficients and memory variables
  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  integer,intent(in) :: iphase

  !---
  !--- local variables
  !---
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxi_m2,dux_dgamma_m2,duz_dxi_m2,duz_dgamma_m2 !lucas, CTD-SEM

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_m2,dux_dzl_m2,duz_dxl_m2,duz_dzl_m2 !lucas, CTD-SEM
 ! real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: dummy_loc_m2 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1_m2,tempx2_m2,tempz1_m2,tempz2_m2 !lucas, CTD-SEM
 ! real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: tempx3
 ! real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: sigma_thetatheta
 ! real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  real(kind=CUSTOM_REAL) :: tempx1l_m2,tempx2l_m2,tempz1l_m2,tempz2l_m2 !lucas, CTD-SEM
  !real(kind=CUSTOM_REAL) :: theta,ct,st
  real(kind=CUSTOM_REAL) :: sigma_xx_m2,sigma_xy_m2,sigma_xz_m2,sigma_zy_m2,sigma_zz_m2,sigma_zx_m2 !lucas, CTD-SEM
 ! real(kind=CUSTOM_REAL) :: sigma_xx_prime,sigma_xz_prime,sigma_zz_prime,sigma_zx_prime
 ! real(kind=CUSTOM_REAL) :: xxi


  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: e1_sum_m2,e11_sum_m2,e13_sum_m2 !lucas, CTD-SEM

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic_m2,lambdal_unrelaxed_elastic_m2, &
    lambdaplus2mu_unrelaxed_elastic_m2,cpl_m2,csl_m2,rhol_m2,lambdalplusmul_unrelaxed_elastic_m2  !lucas, CTD-SEM

  ! for anisotropy
 ! double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25,c22

  ! CPML coefficients and memory variables
 ! real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: accel_elastic_PML
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  !additional variables for CPML implementation in viscoelastic simulation
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl,mu_duz_dzl, &
 !                                                   mu_dux_dzl,mu_duz_dxl
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
 !                                                   mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl,mu_pml_duz_dxl
   double precision :: qkappal,qmul
   double precision :: qkappal_m2,qmul_m2 !lucas, CTD-SEM

  integer :: num_elements,ispec_p

  ! this to avoid a warning at execution time about an undefined variable being used
  ! for the SH component in the case of a P-SV calculation, and vice versa
  !lucas, CTD-SEM
  ! P_SV-case
  sigma_xx_m2 = 0._CUSTOM_REAL
  sigma_xz_m2 = 0._CUSTOM_REAL
  sigma_zz_m2 = 0._CUSTOM_REAL
  sigma_zx_m2 = 0._CUSTOM_REAL
  ! SH-case
  sigma_xy_m2 = 0._CUSTOM_REAL
  sigma_zy_m2 = 0._CUSTOM_REAL

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  ! loop over spectral elements, lucas: to the end of the routine
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! only for elastic spectral elements
    if (.not. ispec_is_elastic(ispec)) cycle
    ! 1.lucas
    if (ATTENUATION_VISCOELASTIC) then
      if (.not. assign_external_model) then
        qkappal = QKappa_attenuation(kmato(ispec))
        qmul = Qmu_attenuation(kmato(ispec))
      else
        qkappal_m2 = maxval(QKappa_attenuationext_m2(:,:,ispec))
        qmul_m2 =  maxval(Qmu_attenuationext_m2(:,:,ispec))
      endif
    endif

    ! 2.lucas
    ! gets local displacement for element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        !lucas1, CTD-SEM
        dummy_loc_m2(1,i,j) = displ_elastic_m2(1,iglob)
        dummy_loc_m2(2,i,j) = displ_elastic_m2(2,iglob)

       !print *, 'displ_elastic_m2(1,10000)=========', displ_elastic_m2(1,10000)
       !print *, 'displ_elastic_m2(2,10000)=========', displ_elastic_m2(2,10000)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
      enddo
    enddo
    !3.lucas
    ! first double loop over GLL points to compute and store gradients
    !call mxm_4comp_singleA(dux_dxi,duz_dxi,dux_dgamma,duz_dgamma,dummy_loc,hprime_xx,hprime_zz) ! lucas: mxm_4comp_singleA(x1,x2,z1,z2,A,B,C), where out: x1,x2,z1,z2.
    !lucas2, CTD-SEM
    call mxm_4comp_singleA(dux_dxi_m2,duz_dxi_m2,dux_dgamma_m2,duz_dgamma_m2,dummy_loc_m2,hprime_xx,hprime_zz) !lucas, CTD-SEM
    
    !4.lucas
    ! AXISYM case overwrites dux_dxi and duz_dxi
   ! if (AXISYM) then
   !   if (is_on_the_axis(ispec)) then
   !     do j = 1,NGLLZ
   !       do i = 1,NGLLX
   !         ! derivative along x and along z
   !         dux_dxi(i,j) = 0._CUSTOM_REAL
   !         duz_dxi(i,j) = 0._CUSTOM_REAL
   !         do k = 1,NGLJ
   !           dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
   !           duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprimeBar_xx(i,k)
   !         enddo
   !       enddo
   !     enddo
   !   endif
   ! endif

    ! 5.lucas
    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        ! lucas3, CTD-SEM
        ! derivatives of displacement
        dux_dxl_m2(i,j) = dux_dxi_m2(i,j)*xixl + dux_dgamma_m2(i,j)*gammaxl
        dux_dzl_m2(i,j) = dux_dxi_m2(i,j)*xizl + dux_dgamma_m2(i,j)*gammazl

        duz_dxl_m2(i,j) = duz_dxi_m2(i,j)*xixl + duz_dgamma_m2(i,j)*gammaxl
        duz_dzl_m2(i,j) = duz_dxi_m2(i,j)*xizl + duz_dgamma_m2(i,j)*gammazl
      enddo
    enddo

    !print *,'dux_dxl_m2(2,2)==',dux_dxl_m2(2,2)
    !print *,'dux_dzl_m2(2,2)==',dux_dzl_m2(2,2)
    !print *,'duz_dxl_m2(2,2)==',duz_dxl_m2(2,2)
    !print *,'duz_dzl_m2(2,2)==',duz_dzl_m2(2,2)

    !6.lucas
    ! AXISYM case overwrite duz_dxl
 !   if (AXISYM) then
 !     if (is_on_the_axis(ispec)) then
 !       ! d_uz/dr=0 on the axis
 !       ! i == 1
 !       do j = 1,NGLLZ
 !         duz_dxl(1,j) = 0._CUSTOM_REAL
 !       enddo
 !     endif
 !   endif
    !7.lucas
 !   if (PML_BOUNDARY_CONDITIONS) then !viscoelastic PML
 !     if (ispec_is_PML(ispec)) then
 !       if (ATTENUATION_VISCOELASTIC) then
 !         if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then

 !          call pml_compute_memory_variables_viscoelastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
 !                                                          kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
 !                                                          mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl, &
 !                                                          mu_pml_duz_dxl,kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl, &
 !                                                          mu_duz_dzl,mu_dux_dzl,mu_duz_dxl)
 !         else
 !           call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
 !                                                    dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
 !                                                    PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
 !                                                    PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
 !         endif
 !       else
 !         call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
 !                                                  dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
 !                                                  PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
 !                                                  PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
 !       endif
 !     endif ! ispec_is_PML
 !   endif ! PML_BOUNDARY_CONDITIONS


    !8.lucas
    ! AXISYM case overwrite duz_dxl and duz_dxl_prime
 !   if (AXISYM) then
 !     if (is_on_the_axis(ispec)) then
 !       ! d_uz/dr=0 on the axis
 !       !i == 1
 !       do j = 1,NGLLZ
 !         duz_dxl(1,j) = 0._CUSTOM_REAL
 !         duz_dxl_prime(1,j) = 0._CUSTOM_REAL
 !       enddo
 !     endif
 !   endif
    !9.lucas
    ! get unrelaxed elastic parameters of current spectral element

    !lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    !mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
    !lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))
    !lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic

    !10.lucas
    ! first double loop to compute gradient
    do j = 1,NGLLZ
      do i = 1,NGLLX
        !--- if external medium, get elastic parameters of current grid point
        ! 10.1 lucas
        ! lucas4, CTD-SEM
        if (assign_external_model) then
          cpl_m2 = vpext_m2(i,j,ispec)
          csl_m2 = vsext_m2(i,j,ispec)
          rhol_m2 = rhoext_m2(i,j,ispec)

          mul_unrelaxed_elastic_m2 = rhol_m2*csl_m2*csl_m2
          lambdal_unrelaxed_elastic_m2 = rhol_m2*cpl_m2*cpl_m2 - 2._CUSTOM_REAL * mul_unrelaxed_elastic_m2
          lambdaplus2mu_unrelaxed_elastic_m2 = lambdal_unrelaxed_elastic_m2 + 2._CUSTOM_REAL * mul_unrelaxed_elastic_m2
          lambdalplusmul_unrelaxed_elastic_m2 = lambdal_unrelaxed_elastic_m2 + mul_unrelaxed_elastic_m2
        endif
         ! print *,'lambdalplusmul_unrelaxed_elastic_m2==',lambdalplusmul_unrelaxed_elastic_m2
        ! 10.2 lucas
        ! compute stress tensor (include attenuation or anisotropy if needed)
       ! if (AXISYM) then
       !   if (is_on_the_axis(ispec)) then
       !     if (i == 1) then
       !       ! First GLJ point
       !       sigma_xx = 0._CUSTOM_REAL
       !       sigma_zz = 0._CUSTOM_REAL
       !       sigma_thetatheta(i,j) = 0._CUSTOM_REAL
       !       xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
       !       r_xiplus1(i,j) = xxi
       !       do k = 1,NGLJ
       !         sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
       !         sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
       !         sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
       !       enddo
       !       sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
       !                 + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + sigma_xx/xxi)
       !       sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
       !                + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + sigma_zz/xxi)
       !       sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
       !       sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
       !                          + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta(i,j)/xxi
       !     else
       !       ! Not first GLJ point
       !       sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
       !                  + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !       sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
       !                  + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !       sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
       !       sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
       !                               + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
       !       r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
       !     endif
       !   else
       !     ! Not on the axis
       !     sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
       !                + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !     sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
       !                + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !     sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
       !     sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
       !                             + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
       !   endif
       ! else ! Not axisym
           !lucas5, CTD-SEM
          if (P_SV) then
            ! P_SV case
            sigma_xx_m2 = lambdaplus2mu_unrelaxed_elastic_m2 * dux_dxl_m2(i,j) + lambdal_unrelaxed_elastic_m2 * duz_dzl_m2(i,j)
            sigma_zz_m2 = lambdaplus2mu_unrelaxed_elastic_m2 * duz_dzl_m2(i,j) + lambdal_unrelaxed_elastic_m2 * dux_dxl_m2(i,j)
            sigma_xz_m2 = mul_unrelaxed_elastic_m2 * (duz_dxl_m2(i,j) + dux_dzl_m2(i,j))
            sigma_zx_m2 = sigma_xz_m2
          else
            ! SH-case
            sigma_xy_m2 = mul_unrelaxed_elastic_m2 * dux_dxl_m2(i,j)
            sigma_zy_m2 = mul_unrelaxed_elastic_m2 * dux_dzl_m2(i,j)
          endif
          !print *,'sigma_zx_m2==',sigma_zx_m2
       ! endif ! AXISYM

        ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
        ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below.

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

        ! 10.3 lucas
        if (ATTENUATION_VISCOELASTIC) then
          ! This routine updates the memory variables and with the current grad(displ)
          ! and gets the attenuation contribution (in e*_sum variables)
          call compute_attenuation_viscoelastic_m2(e1_m2,e11_m2,e13_m2,dux_dxl_m2(i,j),dux_dzl_m2(i,j), &
                                                duz_dxl_m2(i,j),duz_dzl_m2(i,j), &
                                                dux_dxl_old_m2,duz_dzl_old_m2, &
                                                dux_dzl_plus_duz_dxl_old_m2,PML_BOUNDARY_CONDITIONS,i,j,ispec, &
                                                e1_sum_m2,e11_sum_m2,e13_sum_m2) !lucas, CTD-SEM

! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
          sigma_xx_m2 = sigma_xx_m2 + lambdalplusmul_unrelaxed_elastic_m2 * e1_sum_m2 + TWO * mul_unrelaxed_elastic_m2 * e11_sum_m2
          sigma_xz_m2 = sigma_xz_m2 + mul_unrelaxed_elastic_m2 * e13_sum_m2
          sigma_zz_m2 = sigma_zz_m2 + lambdalplusmul_unrelaxed_elastic_m2 * e1_sum_m2 - TWO * mul_unrelaxed_elastic_m2 * e11_sum_m2
          sigma_zx_m2 = sigma_xz_m2
        endif ! ATTENUATION_VISCOELASTIC

        ! 10.4 lucas
      !  if (PML_BOUNDARY_CONDITIONS) then
      !    if (ATTENUATION_VISCOELASTIC) then
      !      if (ispec_is_PML(ispec)) then
      !         if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then
      !           sigma_xx = (lambdalplusmul_unrelaxed_elastic * kappa_pml_dux_dxl(i,j) + mul_unrelaxed_elastic * &
      !                      mu_pml_dux_dxl(i,j)) + (lambdalplusmul_unrelaxed_elastic * kappa_duz_dzl(i,j) - &
      !                      mul_unrelaxed_elastic * mu_duz_dzl(i,j))
      !           sigma_xz = (mul_unrelaxed_elastic * mu_pml_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_dux_dzl(i,j))
      !           sigma_zx = (mul_unrelaxed_elastic * mu_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_pml_dux_dzl(i,j))
      !           sigma_zz = (lambdalplusmul_unrelaxed_elastic * kappa_dux_dxl(i,j) - mul_unrelaxed_elastic * &
      !                      mu_dux_dxl(i,j))+ (lambdalplusmul_unrelaxed_elastic * kappa_pml_duz_dzl(i,j) + &
      !                      mul_unrelaxed_elastic * mu_pml_duz_dzl(i,j))
      !         else
      !           sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * PML_duz_dzl(i,j)
      !           sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * PML_dux_dxl(i,j)
      !           sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
      !           sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
      !         endif
      !      endif
      !    else ! No ATTENUATION
      !      if (ispec_is_PML(ispec) .and. nspec_PML > 0) then
      !        if (ROTATE_PML_ACTIVATE) then
      !          theta = - ROTATE_PML_ANGLE/180._CUSTOM_REAL*PI
      !          if (it == 1) write(*,*) theta,ROTATE_PML_ACTIVATE,cos(theta),sin(theta)
      !          ct = cos(theta)
      !          st = sin(theta)
      !          sigma_xx_prime = lambdaplus2mu_unrelaxed_elastic * (ct**2*dux_dxl(i,j) + ct*st*duz_dxl(i,j) + &
      !                                                              ct*st*dux_dzl(i,j) + st**2*duz_dzl(i,j)) &
      !                           + lambdal_unrelaxed_elastic*(st**2*PML_dux_dxl(i,j) - ct*st*PML_duz_dxl(i,j) - &
      !                                                        ct*st*PML_dux_dzl(i,j) + ct**2*PML_duz_dzl(i,j))
      !
      !          sigma_xz_prime = mul_unrelaxed_elastic * (-ct*st*dux_dxl(i,j) + ct**2*duz_dxl(i,j) - &
      !                                                    st**2*dux_dzl(i,j) + ct*st*duz_dzl(i,j)) &
      !                           + mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) - st**2*PML_duz_dxl(i,j) + &
      !                                                      ct**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j))
      !
      !          sigma_zx_prime = mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) + ct**2*PML_duz_dxl(i,j) - &
      !                                                    st**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j)) &
      !                          + mul_unrelaxed_elastic * (-ct*st*dux_dxl_prime(i,j) - st**2*duz_dxl_prime(i,j) + &
      !                                                    ct**2*dux_dzl_prime(i,j) + ct*st*duz_dzl_prime(i,j))
      !
      !          sigma_zz_prime = lambdaplus2mu_unrelaxed_elastic*(st**2*dux_dxl_prime(i,j) - ct*st*duz_dxl_prime(i,j) - &
      !                                                            ct*st*dux_dzl_prime(i,j) + ct**2*duz_dzl_prime(i,j)) &
      !                           + lambdal_unrelaxed_elastic*(ct**2*PML_dux_dxl(i,j) + ct*st*PML_duz_dxl(i,j) + &
      !                                                        ct*st*PML_dux_dzl(i,j) + st**2*PML_duz_dzl(i,j))
      ! 
      !          sigma_xx = ct**2*sigma_xx_prime - ct*st*sigma_xz_prime - ct*st*sigma_zx_prime + st**2*sigma_zz_prime
      !          sigma_xz = ct*st*sigma_xx_prime + ct**2*sigma_xz_prime - st**2*sigma_zx_prime - ct*st*sigma_zz_prime
      !          sigma_zx = ct*st*sigma_xx_prime - st**2*sigma_xz_prime + ct**2*sigma_zx_prime - ct*st*sigma_zz_prime
      !          sigma_zz = st**2*sigma_xx_prime + ct*st*sigma_xz_prime + ct*st*sigma_zx_prime + ct**2*sigma_zz_prime
      !        else
      !          sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
      !          sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
      !          sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
      !          sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
      !        endif
      !      endif
      !
      !    endif ! ATTENUATION
      !  endif ! PML_BOUNDARY_CONDITION
        ! 10.5 lucas
        ! full anisotropy
      !  if (ispec_is_anisotropic(ispec)) then
      !    if (assign_external_model) then
      !      c11 = c11ext(i,j,ispec)
      !      c13 = c13ext(i,j,ispec)
      !      c15 = c15ext(i,j,ispec)
      !      c33 = c33ext(i,j,ispec)
      !      c35 = c35ext(i,j,ispec)
      !      c55 = c55ext(i,j,ispec)
      !      c12 = c12ext(i,j,ispec)
      !      c23 = c23ext(i,j,ispec)
      !      c25 = c25ext(i,j,ispec)
      !      if (AXISYM) then
      !        c22 = c22ext(i,j,ispec) ! This variable is used for axisym simulations only
      !      endif
      !    else
      !      c11 = anisotropy(1,kmato(ispec))
      !      c13 = anisotropy(2,kmato(ispec))
      !      c15 = anisotropy(3,kmato(ispec))
      !      c33 = anisotropy(4,kmato(ispec))
      !      c35 = anisotropy(5,kmato(ispec))
      !      c55 = anisotropy(6,kmato(ispec))
      !      c12 = anisotropy(7,kmato(ispec))
      !      c23 = anisotropy(8,kmato(ispec))
      !      c25 = anisotropy(9,kmato(ispec))
      !      if (AXISYM) then
      !        c22 = anisotropy(10,kmato(ispec)) ! This variable is used for axisym simulations only
      !      endif
      !    endif

      !    ! implement anisotropy in 2D
      !    if (AXISYM) then
      !      if (is_on_the_axis(ispec)) then
      !        if (i == 1) then
      !          ! first GLJ point, on the axis
      !          sigma_xx = 0._CUSTOM_REAL
      !          sigma_zz = 0._CUSTOM_REAL
      !          sigma_thetatheta(i,j) = 0._CUSTOM_REAL
      !          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
      !          r_xiplus1(i,j) = xxi
      !          do k = 1,NGLJ ! Compute the sum
      !            sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
      !            sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
      !            sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
      !          enddo
      !          sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*sigma_xx/xxi
      !          sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*sigma_zz/xxi
      !          sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !          sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*sigma_thetatheta(i,j)/xxi

      !       else
      !          ! not first GLJ point but not on the axis
      !          sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !          sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !          sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !          sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + &
      !                                  c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !        endif
      !      else
      !        ! axisym but not on the axis
      !        sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !        sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !        sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !        sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !      endif
      !    else
      !      ! not AXISYM
      !      sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c15*(duz_dxl(i,j) + dux_dzl(i,j))
      !      sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c35*(duz_dxl(i,j) + dux_dzl(i,j))
      !      sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !      sigma_zx = sigma_xz
      !    endif
      !  endif
        
        ! 10.6 lucas
        ! weak formulation term based on stress tensor (non-symmetric form)
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        jacobianl = deriv(5,i,j)

        !! ABAB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
        ! tempx1(i,j) = w.J.F_{11}^{ij}
        ! tempz1(i,j) = w.J.F_{21}^{ij}
        ! tempx2(i,j) = w.J.F_{12}^{ij}
        ! tempz2(i,j) = w.J.F_{22}^{ij}
        ! 10.7 lucas
      !  if (AXISYM) then
      !    if (is_on_the_axis(ispec)) then
      !      ! This is normal, we always add a contribution depending on the value on the axis
      !      ! i.e. we purposely sum something at point (i,j) with something at point (1,j)
      !      tempx3(i,j) = jacobian(1,j,ispec) * sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)
      
      !      ! not first GLJ point
      !      if (abs(coord(1,ibool(i,j,ispec))) > TINYVAL) then
      !        if (i == 1) then
      !          write(*,*) "Element number:",ispec
      !          call stop_the_code("Error: an axial element is rotated. The code should have been stopped before. Check that your &
      !           &coordinates are greater than TINYVAL. Maybe you should also have a look to &
      !           &doc/problematic_case_that_we_exclude_for_axisymmetric.pdf")
      !        endif
      !        tempx3(i,j) = tempx3(i,j) + wxglj(i) * jacobianl &
      !                      * sigma_thetatheta(i,j)/(xiglj(i)+ONE) ! this goes to accel_x
      !      endif

      !     tempx2(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
      !      tempz2(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
      !      tempx1(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
      !      tempz1(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
      !    else
      !      ! axisym but not on the axis
      !      tempx2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
      !      tempz2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
      !      tempx1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
      !      tempz1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
      !      tempx3(i,j) = wxgll(i) * jacobianl * sigma_thetatheta(i,j) ! this goes to accel_x
      !    endif
      !  else
          !lucas6, CTD-SEM
          ! default (not axisym case)

   

          if (P_SV) then
            ! P_SV case
            tempx1_m2(i,j) = jacobianl * (sigma_xx_m2*xixl + sigma_zx_m2*xizl) ! this goes to accel_x
            tempz1_m2(i,j) = jacobianl * (sigma_xz_m2*xixl + sigma_zz_m2*xizl) ! this goes to accel_z

            tempx2_m2(i,j) = jacobianl * (sigma_xx_m2*gammaxl + sigma_zx_m2*gammazl) ! this goes to accel_x
            tempz2_m2(i,j) = jacobianl * (sigma_xz_m2*gammaxl + sigma_zz_m2*gammazl) ! this goes to accel_z

          !print *,'Lucas test here: '
          !print *,'jacobianl = ', jacobianl
          !print *,'xixl = ', xixl
          !print *,'xizl = ', xizl
          !print *,'gammaxl = ', gammaxl
          !print *,'gammazl = ', gammazl

          !print *,'tempx1_m2(2,2) = ', tempx1_m2(2,2)
          !print *,'tempz1_m2(2,2) = ', tempz1_m2(2,2)
          !print *,'tempx2_m2(2,2) = ', tempx2_m2(2,2)
          !print *,'tempz2_m2(2,2) = ', tempz2_m2(2,2)
       
 
          else
            ! SH-case
            tempx1_m2(i,j) = jacobianl * (sigma_xy_m2*xixl + sigma_zy_m2*xizl) ! this goes to accel_x
            tempx2_m2(i,j) = jacobianl * (sigma_xy_m2*gammaxl + sigma_zy_m2*gammazl) ! this goes to accel_x
            tempz1_m2(i,j) = 0._CUSTOM_REAL
            tempz2_m2(i,j) = 0._CUSTOM_REAL
          endif
      !  endif ! AXISYM
      enddo
    enddo  ! end of the loops on the collocation points i,j
    !11.lucas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update the displacement memory variable
!    if (PML_BOUNDARY_CONDITIONS) then
!      ! calculates contribution from each C-PML element to update acceleration
!      call pml_compute_accel_contribution_elastic(ispec,nglob, &
!                                                   dummy_loc,displ_elastic_old,veloc_elastic, &
!                                                   accel_elastic_PML,r_xiplus1)
!    endif
    !12.lucas
    !
    ! second double-loop over GLL to compute all the terms
    ! 12.1 lucas
  !  if (AXISYM) then
  !    if (is_on_the_axis(ispec)) then

  !      do j = 1,NGLLZ
  !        do i = 1,NGLLX
  !          iglob = ibool(i,j,ispec)
  !          if (.not. iglob_is_forced(iglob)) then
  !            ! along x direction and z direction
  !            ! and assemble the contributions
  !            ! we can merge the two loops because NGLLX == NGLLZ

  !            ! assembles the contributions
  !            tempx1l = 0._CUSTOM_REAL
  !            tempx2l = 0._CUSTOM_REAL
  !            tempz1l = 0._CUSTOM_REAL
  !            tempz2l = 0._CUSTOM_REAL
  !            do k = 1,NGLLX
  !              tempx1l = tempx1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
  !              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
  !              tempz1l = tempz1l + tempz1(k,j) * hprimeBarwglj_xx(k,i)
  !              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
  !            enddo
  !            ! sums contributions from each element to the global values
  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxglj(i) * tempx2l)
  !            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxglj(i) * tempz2l)

  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
  !          endif
  !        enddo
  !      enddo

  !    else
  !      ! Axisym but not on the axis
  !      do j = 1,NGLLZ
  !        do i = 1,NGLLX
  !          iglob = ibool(i,j,ispec)
  !          if (.not. iglob_is_forced(iglob)) then
  !            ! assembles the contributions
  !            tempx1l = 0._CUSTOM_REAL
  !            tempx2l = 0._CUSTOM_REAL
  !            tempz1l = 0._CUSTOM_REAL
  !            tempz2l = 0._CUSTOM_REAL
  !            do k = 1,NGLLX
  !              tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
  !              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
  !              tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
  !              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
  !            enddo
  !            ! sums contributions from each element to the global values
  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
  !            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)

  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
  !          endif
  !        enddo
  !      enddo
  !    endif
  !  else ! 12.2 lucas
      !lucas7, CTD-SEM
      ! default (not AXISYM case)
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (.not. iglob_is_forced(iglob)) then
            ! assembles the contributions
            tempx1l_m2 = 0._CUSTOM_REAL
            tempx2l_m2 = 0._CUSTOM_REAL
            tempz1l_m2 = 0._CUSTOM_REAL
            tempz2l_m2 = 0._CUSTOM_REAL
            do k = 1,NGLLX
              tempx1l_m2 = tempx1l_m2 + tempx1_m2(k,j) * hprimewgll_xx(k,i)
              tempx2l_m2 = tempx2l_m2 + tempx2_m2(i,k) * hprimewgll_zz(k,j)
              tempz1l_m2 = tempz1l_m2 + tempz1_m2(k,j) * hprimewgll_xx(k,i)
              tempz2l_m2 = tempz2l_m2 + tempz2_m2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            accel_elastic_m2(1,iglob) = accel_elastic_m2(1,iglob) - (wzgll(j) * tempx1l_m2 + wxgll(i) * tempx2l_m2)
            accel_elastic_m2(2,iglob) = accel_elastic_m2(2,iglob) - (wzgll(j) * tempz1l_m2 + wxgll(i) * tempz2l_m2)
          endif
        enddo
      enddo
  !  endif ! AXISYM

    !13.lucas
    !lucas8, CTD-SEM
    ! adds PML_BOUNDARY_CONDITIONS contribution
 !   if (PML_BOUNDARY_CONDITIONS) then
 !     if (ispec_is_PML(ispec)) then
 !       do j = 1,NGLLZ
 !         do i = 1,NGLLX
 !           iglob = ibool(i,j,ispec)
 !           if (.not. iglob_is_forced(iglob)) then
 !             accel_elastic(1,iglob) = accel_elastic(1,iglob) - accel_elastic_PML(1,i,j)
 !             accel_elastic(2,iglob) = accel_elastic(2,iglob) - accel_elastic_PML(2,i,j)
 !           endif
 !         enddo
 !       enddo
 !     endif
 !   endif

  enddo ! end of loop over all spectral elements

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_4comp_singleA(x1,x2,z1,z2,A,B,C)

! matrix x matrix multiplication, merging 4 loops for x1,x2 = A^t B^t and z1,z2 = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          duz_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!          duz_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprime_xx(i,k)
!            duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + dummy_loc(1,i,k)*hprime_zz(j,k)
!            duz_dgamma(i,j) = duz_dgamma(i,j) + dummy_loc(2,i,k)*hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NDIM,NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x1,x2,z1,z2

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: A
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x1(i,j) = A(1,1,j) * B(i,1) + A(1,2,j) * B(i,2) + A(1,3,j) * B(i,3) + A(1,4,j) * B(i,4) + A(1,5,j) * B(i,5)
        x2(i,j) = A(2,1,j) * B(i,1) + A(2,2,j) * B(i,2) + A(2,3,j) * B(i,3) + A(2,4,j) * B(i,4) + A(2,5,j) * B(i,5)

        z1(i,j) = A(1,i,1) * C(j,1) + A(1,i,2) * C(j,2) + A(1,i,3) * C(j,3) + A(1,i,4) * C(j,4) + A(1,i,5) * C(j,5)
        z2(i,j) = A(2,i,1) * C(j,1) + A(2,i,2) * C(j,2) + A(2,i,3) * C(j,3) + A(2,i,4) * C(j,4) + A(2,i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x1(i,j) = 0._CUSTOM_REAL
        x2(i,j) = 0._CUSTOM_REAL
        z1(i,j) = 0._CUSTOM_REAL
        z2(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x1(i,j) = x1(i,j) + A(1,k,j) * B(i,k)
          x2(i,j) = x2(i,j) + A(2,k,j) * B(i,k)

          z1(i,j) = z1(i,j) + A(1,i,k) * C(j,k)
          z2(i,j) = z2(i,j) + A(2,i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_4comp_singleA

  end subroutine compute_forces_viscoelastic_m2 ! lucas, CTD-SEM


!! ----------------------------------------------------------m1 -------------------------------------------------------------
!! Yujiang Lucas Xie

  subroutine compute_forces_viscoelastic_m1(accel_elastic_m1,veloc_elastic_m1,displ_elastic_m1,displ_elastic_old_m1, &
                                         dux_dxl_old_m1,duz_dzl_old_m1,dux_dzl_plus_duz_dxl_old_m1, &
                                         PML_BOUNDARY_CONDITIONS,e1_m1,e11_m1,e13_m1,iphase) ! lucas, CTD-SEM
                                                                                     

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    ONE,TWO,PI,TINYVAL,FOUR_THIRDS,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: nglob,assign_external_model,P_SV, &
                         ATTENUATION_VISCOELASTIC,nspec_ATT_el,N_SLS, &
                         ibool,kmato,ispec_is_elastic, &
                         xix,xiz,gammax,gammaz, &
                         jacobian,vpext,vsext,rhoext, & ! lucas, CTD-SEM
                         QKappa_attenuation,Qmu_attenuation,QKappa_attenuationext,Qmu_attenuationext, & ! lucas, CTD-SEM
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         iglob_is_forced

  ! overlapping communication
  use specfem_par, only: nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  ! AXISYM
 ! use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  ! PML arrays
 ! use specfem_par, only: nspec_PML,ispec_is_PML,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_m1,veloc_elastic_m1 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic_m1 !lucas, CTD-SEM

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_old_m1 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: dux_dxl_old_m1,duz_dzl_old_m1 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: dux_dzl_plus_duz_dxl_old_m1 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el,N_SLS),intent(inout) :: e1_m1,e11_m1,e13_m1 !lucas, CTD-SEM

  ! CPML coefficients and memory variables
  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  integer,intent(in) :: iphase

  !---
  !--- local variables
  !---
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxi_m1,dux_dgamma_m1,duz_dxi_m1,duz_dgamma_m1 !lucas, CTD-SEM

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_m1,dux_dzl_m1,duz_dxl_m1,duz_dzl_m1 !lucas, CTD-SEM
 ! real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: dummy_loc_m1 !lucas, CTD-SEM
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1_m1,tempx2_m1,tempz1_m1,tempz2_m1 !lucas, CTD-SEM
 ! real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: tempx3
 ! real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: sigma_thetatheta
 ! real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  real(kind=CUSTOM_REAL) :: tempx1l_m1,tempx2l_m1,tempz1l_m1,tempz2l_m1 !lucas, CTD-SEM
  !real(kind=CUSTOM_REAL) :: theta,ct,st
  real(kind=CUSTOM_REAL) :: sigma_xx_m1,sigma_xy_m1,sigma_xz_m1,sigma_zy_m1,sigma_zz_m1,sigma_zx_m1 !lucas, CTD-SEM
 ! real(kind=CUSTOM_REAL) :: sigma_xx_prime,sigma_xz_prime,sigma_zz_prime,sigma_zx_prime
 ! real(kind=CUSTOM_REAL) :: xxi


  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: e1_sum_m1,e11_sum_m1,e13_sum_m1 !lucas, CTD-SEM

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic_m1,lambdal_unrelaxed_elastic_m1, &
    lambdaplus2mu_unrelaxed_elastic_m1,cpl_m1,csl_m1,rhol_m1,lambdalplusmul_unrelaxed_elastic_m1  !lucas, CTD-SEM

  ! for anisotropy
 ! double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25,c22

  ! CPML coefficients and memory variables
 ! real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: accel_elastic_PML
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  !additional variables for CPML implementation in viscoelastic simulation
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl,mu_duz_dzl, &
 !                                                   mu_dux_dzl,mu_duz_dxl
 ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
 !                                                   mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl,mu_pml_duz_dxl
  double precision :: qkappal,qmul

  integer :: num_elements,ispec_p

  ! this to avoid a warning at execution time about an undefined variable being used
  ! for the SH component in the case of a P-SV calculation, and vice versa
  !lucas, CTD-SEM
  ! P_SV-case
  sigma_xx_m1 = 0._CUSTOM_REAL
  sigma_xz_m1 = 0._CUSTOM_REAL
  sigma_zz_m1 = 0._CUSTOM_REAL
  sigma_zx_m1 = 0._CUSTOM_REAL
  ! SH-case
  sigma_xy_m1 = 0._CUSTOM_REAL
  sigma_zy_m1 = 0._CUSTOM_REAL

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  ! loop over spectral elements, lucas: to the end of the routine
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! only for elastic spectral elements
    if (.not. ispec_is_elastic(ispec)) cycle
    ! 1.lucas
    if (ATTENUATION_VISCOELASTIC) then
      if (.not. assign_external_model) then
        qkappal = QKappa_attenuation(kmato(ispec))
        qmul = Qmu_attenuation(kmato(ispec))
      else
        qkappal = maxval(QKappa_attenuationext(:,:,ispec))
        qmul =  maxval(Qmu_attenuationext(:,:,ispec))
      endif
    endif

    ! 2.lucas
    ! gets local displacement for element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        !lucas1, CTD-SEM
        dummy_loc_m1(1,i,j) = displ_elastic_m1(1,iglob)
        dummy_loc_m1(2,i,j) = displ_elastic_m1(2,iglob)

       !print *, 'displ_elastic_m1(1,10000)=========', displ_elastic_m1(1,10000)
       !print *, 'displ_elastic_m1(2,10000)=========', displ_elastic_m1(2,10000)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
      enddo
    enddo
    !3.lucas
    ! first double loop over GLL points to compute and store gradients
    !call mxm_4comp_singleA(dux_dxi,duz_dxi,dux_dgamma,duz_dgamma,dummy_loc,hprime_xx,hprime_zz) ! lucas: mxm_4comp_singleA(x1,x2,z1,z2,A,B,C), where out: x1,x2,z1,z2.
    !lucas2, CTD-SEM
    call mxm_4comp_singleA(dux_dxi_m1,duz_dxi_m1,dux_dgamma_m1,duz_dgamma_m1,dummy_loc_m1,hprime_xx,hprime_zz) !lucas, CTD-SEM
    
    !4.lucas
    ! AXISYM case overwrites dux_dxi and duz_dxi
   ! if (AXISYM) then
   !   if (is_on_the_axis(ispec)) then
   !     do j = 1,NGLLZ
   !       do i = 1,NGLLX
   !         ! derivative along x and along z
   !         dux_dxi(i,j) = 0._CUSTOM_REAL
   !         duz_dxi(i,j) = 0._CUSTOM_REAL
   !         do k = 1,NGLJ
   !           dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
   !           duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprimeBar_xx(i,k)
   !         enddo
   !       enddo
   !     enddo
   !   endif
   ! endif

    ! 5.lucas
    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        ! lucas3, CTD-SEM
        ! derivatives of displacement
        dux_dxl_m1(i,j) = dux_dxi_m1(i,j)*xixl + dux_dgamma_m1(i,j)*gammaxl
        dux_dzl_m1(i,j) = dux_dxi_m1(i,j)*xizl + dux_dgamma_m1(i,j)*gammazl

        duz_dxl_m1(i,j) = duz_dxi_m1(i,j)*xixl + duz_dgamma_m1(i,j)*gammaxl
        duz_dzl_m1(i,j) = duz_dxi_m1(i,j)*xizl + duz_dgamma_m1(i,j)*gammazl
      enddo
    enddo

    !print *,'dux_dxl_m1(2,2)==',dux_dxl_m1(2,2)
    !print *,'dux_dzl_m1(2,2)==',dux_dzl_m1(2,2)
    !print *,'duz_dxl_m1(2,2)==',duz_dxl_m1(2,2)
    !print *,'duz_dzl_m1(2,2)==',duz_dzl_m1(2,2)

    !6.lucas
    ! AXISYM case overwrite duz_dxl
 !   if (AXISYM) then
 !     if (is_on_the_axis(ispec)) then
 !       ! d_uz/dr=0 on the axis
 !       ! i == 1
 !       do j = 1,NGLLZ
 !         duz_dxl(1,j) = 0._CUSTOM_REAL
 !       enddo
 !     endif
 !   endif

    !7.lucas
 !   if (PML_BOUNDARY_CONDITIONS) then !viscoelastic PML
 !     if (ispec_is_PML(ispec)) then
 !       if (ATTENUATION_VISCOELASTIC) then
 !         if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then

 !          call pml_compute_memory_variables_viscoelastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
 !                                                          kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
 !                                                          mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl, &
 !                                                          mu_pml_duz_dxl,kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl, &
 !                                                          mu_duz_dzl,mu_dux_dzl,mu_duz_dxl)
 !         else
 !           call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
 !                                                    dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
 !                                                    PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
 !                                                    PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
 !         endif
 !       else
 !         call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
 !                                                  dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
 !                                                  PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
 !                                                  PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
 !       endif
 !     endif ! ispec_is_PML
 !   endif ! PML_BOUNDARY_CONDITIONS


    !8.lucas
    ! AXISYM case overwrite duz_dxl and duz_dxl_prime
 !   if (AXISYM) then
 !     if (is_on_the_axis(ispec)) then
 !       ! d_uz/dr=0 on the axis
 !       !i == 1
 !       do j = 1,NGLLZ
 !         duz_dxl(1,j) = 0._CUSTOM_REAL
 !         duz_dxl_prime(1,j) = 0._CUSTOM_REAL
 !       enddo
 !     endif
 !   endif
    !9.lucas
    ! get unrelaxed elastic parameters of current spectral element

    !lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    !mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
    !lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))
    !lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic

    !10.lucas
    ! first double loop to compute gradient
    do j = 1,NGLLZ
      do i = 1,NGLLX
        !--- if external medium, get elastic parameters of current grid point
        ! 10.1 lucas
        ! lucas4, CTD-SEM
        if (assign_external_model) then
          cpl_m1 = vpext(i,j,ispec) !lucas, vpext for m1
          csl_m1 = vsext(i,j,ispec)
          rhol_m1 = rhoext(i,j,ispec)

          mul_unrelaxed_elastic_m1 = rhol_m1*csl_m1*csl_m1
          lambdal_unrelaxed_elastic_m1 = rhol_m1*cpl_m1*cpl_m1 - 2._CUSTOM_REAL * mul_unrelaxed_elastic_m1
          lambdaplus2mu_unrelaxed_elastic_m1 = lambdal_unrelaxed_elastic_m1 + 2._CUSTOM_REAL * mul_unrelaxed_elastic_m1
          lambdalplusmul_unrelaxed_elastic_m1 = lambdal_unrelaxed_elastic_m1 + mul_unrelaxed_elastic_m1
        endif
         ! print *,'lambdalplusmul_unrelaxed_elastic_m1==',lambdalplusmul_unrelaxed_elastic_m1
        ! 10.2 lucas
        ! compute stress tensor (include attenuation or anisotropy if needed)
       ! if (AXISYM) then
       !   if (is_on_the_axis(ispec)) then
       !     if (i == 1) then
       !       ! First GLJ point
       !       sigma_xx = 0._CUSTOM_REAL
       !       sigma_zz = 0._CUSTOM_REAL
       !       sigma_thetatheta(i,j) = 0._CUSTOM_REAL
       !       xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
       !       r_xiplus1(i,j) = xxi
       !       do k = 1,NGLJ
       !         sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
       !         sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
       !         sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
       !       enddo
       !       sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
       !                 + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + sigma_xx/xxi)
       !       sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
       !                + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + sigma_zz/xxi)
       !       sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
       !       sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
       !                          + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta(i,j)/xxi
       !     else
       !       ! Not first GLJ point
       !       sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
       !                  + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !       sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
       !                  + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !       sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
       !       sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
       !                               + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
       !       r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
       !     endif
       !   else
       !     ! Not on the axis
       !     sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
       !                + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !     sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
       !                + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
       !     sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
       !     sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
       !                             + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
       !   endif
       ! else ! Not axisym
           !lucas5, CTD-SEM
          if (P_SV) then
            ! P_SV case
            sigma_xx_m1 = lambdaplus2mu_unrelaxed_elastic_m1 * dux_dxl_m1(i,j) + lambdal_unrelaxed_elastic_m1 * duz_dzl_m1(i,j)
            sigma_zz_m1 = lambdaplus2mu_unrelaxed_elastic_m1 * duz_dzl_m1(i,j) + lambdal_unrelaxed_elastic_m1 * dux_dxl_m1(i,j)
            sigma_xz_m1 = mul_unrelaxed_elastic_m1 * (duz_dxl_m1(i,j) + dux_dzl_m1(i,j))
            sigma_zx_m1 = sigma_xz_m1
          else
            ! SH-case
            sigma_xy_m1 = mul_unrelaxed_elastic_m1 * dux_dxl_m1(i,j)
            sigma_zy_m1 = mul_unrelaxed_elastic_m1 * dux_dzl_m1(i,j)
          endif
          !print *,'sigma_zx_m1==',sigma_zx_m1
       ! endif ! AXISYM

        ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
        ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below.

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

        ! 10.3 lucas
        if (ATTENUATION_VISCOELASTIC) then
          ! This routine updates the memory variables and with the current grad(displ)
          ! and gets the attenuation contribution (in e*_sum variables)
          call compute_attenuation_viscoelastic_m1(e1_m1,e11_m1,e13_m1,dux_dxl_m1(i,j),dux_dzl_m1(i,j), & 
                                                duz_dxl_m1(i,j),duz_dzl_m1(i,j), &
                                                dux_dxl_old_m1,duz_dzl_old_m1, &
                                                dux_dzl_plus_duz_dxl_old_m1,PML_BOUNDARY_CONDITIONS,i,j,ispec, &
                                                e1_sum_m1,e11_sum_m1,e13_sum_m1)

! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
          sigma_xx_m1 = sigma_xx_m1 + lambdalplusmul_unrelaxed_elastic_m1 * e1_sum_m1 + TWO * mul_unrelaxed_elastic_m1 * e11_sum_m1
          sigma_xz_m1 = sigma_xz_m1 + mul_unrelaxed_elastic_m1 * e13_sum_m1
          sigma_zz_m1 = sigma_zz_m1 + lambdalplusmul_unrelaxed_elastic_m1 * e1_sum_m1 - TWO * mul_unrelaxed_elastic_m1 * e11_sum_m1
          sigma_zx_m1 = sigma_xz_m1
        endif ! ATTENUATION_VISCOELASTIC

        ! 10.4 lucas
      !  if (PML_BOUNDARY_CONDITIONS) then
      !    if (ATTENUATION_VISCOELASTIC) then
      !      if (ispec_is_PML(ispec)) then
      !         if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then
      !           sigma_xx = (lambdalplusmul_unrelaxed_elastic * kappa_pml_dux_dxl(i,j) + mul_unrelaxed_elastic * &
      !                      mu_pml_dux_dxl(i,j)) + (lambdalplusmul_unrelaxed_elastic * kappa_duz_dzl(i,j) - &
      !                      mul_unrelaxed_elastic * mu_duz_dzl(i,j))
      !           sigma_xz = (mul_unrelaxed_elastic * mu_pml_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_dux_dzl(i,j))
      !           sigma_zx = (mul_unrelaxed_elastic * mu_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_pml_dux_dzl(i,j))
      !           sigma_zz = (lambdalplusmul_unrelaxed_elastic * kappa_dux_dxl(i,j) - mul_unrelaxed_elastic * &
      !                      mu_dux_dxl(i,j))+ (lambdalplusmul_unrelaxed_elastic * kappa_pml_duz_dzl(i,j) + &
      !                      mul_unrelaxed_elastic * mu_pml_duz_dzl(i,j))
      !         else
      !           sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * PML_duz_dzl(i,j)
      !           sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * PML_dux_dxl(i,j)
      !           sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
      !           sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
      !         endif
      !      endif
      !    else ! No ATTENUATION
      !      if (ispec_is_PML(ispec) .and. nspec_PML > 0) then
      !        if (ROTATE_PML_ACTIVATE) then
      !          theta = - ROTATE_PML_ANGLE/180._CUSTOM_REAL*PI
      !          if (it == 1) write(*,*) theta,ROTATE_PML_ACTIVATE,cos(theta),sin(theta)
      !          ct = cos(theta)
      !          st = sin(theta)
      !          sigma_xx_prime = lambdaplus2mu_unrelaxed_elastic * (ct**2*dux_dxl(i,j) + ct*st*duz_dxl(i,j) + &
      !                                                              ct*st*dux_dzl(i,j) + st**2*duz_dzl(i,j)) &
      !                           + lambdal_unrelaxed_elastic*(st**2*PML_dux_dxl(i,j) - ct*st*PML_duz_dxl(i,j) - &
      !                                                        ct*st*PML_dux_dzl(i,j) + ct**2*PML_duz_dzl(i,j))
      !
      !          sigma_xz_prime = mul_unrelaxed_elastic * (-ct*st*dux_dxl(i,j) + ct**2*duz_dxl(i,j) - &
      !                                                    st**2*dux_dzl(i,j) + ct*st*duz_dzl(i,j)) &
      !                           + mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) - st**2*PML_duz_dxl(i,j) + &
      !                                                      ct**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j))
      !
      !          sigma_zx_prime = mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) + ct**2*PML_duz_dxl(i,j) - &
      !                                                    st**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j)) &
      !                          + mul_unrelaxed_elastic * (-ct*st*dux_dxl_prime(i,j) - st**2*duz_dxl_prime(i,j) + &
      !                                                    ct**2*dux_dzl_prime(i,j) + ct*st*duz_dzl_prime(i,j))
      !
      !          sigma_zz_prime = lambdaplus2mu_unrelaxed_elastic*(st**2*dux_dxl_prime(i,j) - ct*st*duz_dxl_prime(i,j) - &
      !                                                            ct*st*dux_dzl_prime(i,j) + ct**2*duz_dzl_prime(i,j)) &
      !                           + lambdal_unrelaxed_elastic*(ct**2*PML_dux_dxl(i,j) + ct*st*PML_duz_dxl(i,j) + &
      !                                                        ct*st*PML_dux_dzl(i,j) + st**2*PML_duz_dzl(i,j))
      ! 
      !          sigma_xx = ct**2*sigma_xx_prime - ct*st*sigma_xz_prime - ct*st*sigma_zx_prime + st**2*sigma_zz_prime
      !          sigma_xz = ct*st*sigma_xx_prime + ct**2*sigma_xz_prime - st**2*sigma_zx_prime - ct*st*sigma_zz_prime
      !          sigma_zx = ct*st*sigma_xx_prime - st**2*sigma_xz_prime + ct**2*sigma_zx_prime - ct*st*sigma_zz_prime
      !          sigma_zz = st**2*sigma_xx_prime + ct*st*sigma_xz_prime + ct*st*sigma_zx_prime + ct**2*sigma_zz_prime
      !        else
      !          sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
      !          sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
      !          sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
      !          sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
      !        endif
      !      endif
      !
      !    endif ! ATTENUATION
      !  endif ! PML_BOUNDARY_CONDITION
        ! 10.5 lucas
        ! full anisotropy
      !  if (ispec_is_anisotropic(ispec)) then
      !    if (assign_external_model) then
      !      c11 = c11ext(i,j,ispec)
      !      c13 = c13ext(i,j,ispec)
      !      c15 = c15ext(i,j,ispec)
      !      c33 = c33ext(i,j,ispec)
      !      c35 = c35ext(i,j,ispec)
      !      c55 = c55ext(i,j,ispec)
      !      c12 = c12ext(i,j,ispec)
      !      c23 = c23ext(i,j,ispec)
      !      c25 = c25ext(i,j,ispec)
      !      if (AXISYM) then
      !        c22 = c22ext(i,j,ispec) ! This variable is used for axisym simulations only
      !      endif
      !    else
      !      c11 = anisotropy(1,kmato(ispec))
      !      c13 = anisotropy(2,kmato(ispec))
      !      c15 = anisotropy(3,kmato(ispec))
      !      c33 = anisotropy(4,kmato(ispec))
      !      c35 = anisotropy(5,kmato(ispec))
      !      c55 = anisotropy(6,kmato(ispec))
      !      c12 = anisotropy(7,kmato(ispec))
      !      c23 = anisotropy(8,kmato(ispec))
      !      c25 = anisotropy(9,kmato(ispec))
      !      if (AXISYM) then
      !        c22 = anisotropy(10,kmato(ispec)) ! This variable is used for axisym simulations only
      !      endif
      !    endif

      !    ! implement anisotropy in 2D
      !    if (AXISYM) then
      !      if (is_on_the_axis(ispec)) then
      !        if (i == 1) then
      !          ! first GLJ point, on the axis
      !          sigma_xx = 0._CUSTOM_REAL
      !          sigma_zz = 0._CUSTOM_REAL
      !          sigma_thetatheta(i,j) = 0._CUSTOM_REAL
      !          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
      !          r_xiplus1(i,j) = xxi
      !          do k = 1,NGLJ ! Compute the sum
      !            sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
      !            sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
      !            sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
      !          enddo
      !          sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*sigma_xx/xxi
      !          sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*sigma_zz/xxi
      !          sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !          sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*sigma_thetatheta(i,j)/xxi

      !       else
      !          ! not first GLJ point but not on the axis
      !          sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !          sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !          sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !          sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + &
      !                                  c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !        endif
      !      else
      !        ! axisym but not on the axis
      !        sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !        sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !        sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !        sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
      !      endif
      !    else
      !      ! not AXISYM
      !      sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c15*(duz_dxl(i,j) + dux_dzl(i,j))
      !      sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c35*(duz_dxl(i,j) + dux_dzl(i,j))
      !      sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
      !      sigma_zx = sigma_xz
      !    endif
      !  endif
        
        ! 10.6 lucas
        ! weak formulation term based on stress tensor (non-symmetric form)
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        jacobianl = deriv(5,i,j)

        !! ABAB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
        ! tempx1(i,j) = w.J.F_{11}^{ij}
        ! tempz1(i,j) = w.J.F_{21}^{ij}
        ! tempx2(i,j) = w.J.F_{12}^{ij}
        ! tempz2(i,j) = w.J.F_{22}^{ij}
        ! 10.7 lucas
      !  if (AXISYM) then
      !    if (is_on_the_axis(ispec)) then
      !      ! This is normal, we always add a contribution depending on the value on the axis
      !      ! i.e. we purposely sum something at point (i,j) with something at point (1,j)
      !      tempx3(i,j) = jacobian(1,j,ispec) * sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)
      
      !      ! not first GLJ point
      !      if (abs(coord(1,ibool(i,j,ispec))) > TINYVAL) then
      !        if (i == 1) then
      !          write(*,*) "Element number:",ispec
      !          call stop_the_code("Error: an axial element is rotated. The code should have been stopped before. Check that your &
      !           &coordinates are greater than TINYVAL. Maybe you should also have a look to &
      !           &doc/problematic_case_that_we_exclude_for_axisymmetric.pdf")
      !        endif
      !        tempx3(i,j) = tempx3(i,j) + wxglj(i) * jacobianl &
      !                      * sigma_thetatheta(i,j)/(xiglj(i)+ONE) ! this goes to accel_x
      !      endif

      !     tempx2(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
      !      tempz2(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
      !      tempx1(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
      !      tempz1(i,j) = r_xiplus1(i,j) * jacobianl &
      !                    * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
      !    else
      !      ! axisym but not on the axis
      !      tempx2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
      !      tempz2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
      !      tempx1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
      !      tempz1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
      !                    *(sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
      !      tempx3(i,j) = wxgll(i) * jacobianl * sigma_thetatheta(i,j) ! this goes to accel_x
      !    endif
      !  else
          !lucas6, CTD-SEM
          ! default (not axisym case)

   

          if (P_SV) then
            ! P_SV case
            tempx1_m1(i,j) = jacobianl * (sigma_xx_m1*xixl + sigma_zx_m1*xizl) ! this goes to accel_x
            tempz1_m1(i,j) = jacobianl * (sigma_xz_m1*xixl + sigma_zz_m1*xizl) ! this goes to accel_z

            tempx2_m1(i,j) = jacobianl * (sigma_xx_m1*gammaxl + sigma_zx_m1*gammazl) ! this goes to accel_x
            tempz2_m1(i,j) = jacobianl * (sigma_xz_m1*gammaxl + sigma_zz_m1*gammazl) ! this goes to accel_z

          !print *,'Lucas test here: '
          !print *,'jacobianl = ', jacobianl
          !print *,'xixl = ', xixl
          !print *,'xizl = ', xizl
          !print *,'gammaxl = ', gammaxl
          !print *,'gammazl = ', gammazl

          !print *,'tempx1_m1(2,2) = ', tempx1_m1(2,2)
          !print *,'tempz1_m1(2,2) = ', tempz1_m1(2,2)
          !print *,'tempx2_m1(2,2) = ', tempx2_m1(2,2)
          !print *,'tempz2_m1(2,2) = ', tempz2_m1(2,2)
       
 
          else
            ! SH-case
            tempx1_m1(i,j) = jacobianl * (sigma_xy_m1*xixl + sigma_zy_m1*xizl) ! this goes to accel_x
            tempx2_m1(i,j) = jacobianl * (sigma_xy_m1*gammaxl + sigma_zy_m1*gammazl) ! this goes to accel_x
            tempz1_m1(i,j) = 0._CUSTOM_REAL
            tempz2_m1(i,j) = 0._CUSTOM_REAL
          endif
      !  endif ! AXISYM
      enddo
    enddo  ! end of the loops on the collocation points i,j
    !11.lucas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update the displacement memory variable
!    if (PML_BOUNDARY_CONDITIONS) then
!      ! calculates contribution from each C-PML element to update acceleration
!      call pml_compute_accel_contribution_elastic(ispec,nglob, &
!                                                   dummy_loc,displ_elastic_old,veloc_elastic, &
!                                                   accel_elastic_PML,r_xiplus1)
!    endif
    !12.lucas
    !
    ! second double-loop over GLL to compute all the terms
    ! 12.1 lucas
  !  if (AXISYM) then
  !    if (is_on_the_axis(ispec)) then

  !      do j = 1,NGLLZ
  !        do i = 1,NGLLX
  !          iglob = ibool(i,j,ispec)
  !          if (.not. iglob_is_forced(iglob)) then
  !            ! along x direction and z direction
  !            ! and assemble the contributions
  !            ! we can merge the two loops because NGLLX == NGLLZ

  !            ! assembles the contributions
  !            tempx1l = 0._CUSTOM_REAL
  !            tempx2l = 0._CUSTOM_REAL
  !            tempz1l = 0._CUSTOM_REAL
  !            tempz2l = 0._CUSTOM_REAL
  !            do k = 1,NGLLX
  !              tempx1l = tempx1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
  !              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
  !              tempz1l = tempz1l + tempz1(k,j) * hprimeBarwglj_xx(k,i)
  !              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
  !            enddo
  !            ! sums contributions from each element to the global values
  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxglj(i) * tempx2l)
  !            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxglj(i) * tempz2l)

  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
  !          endif
  !        enddo
  !      enddo

  !    else
  !      ! Axisym but not on the axis
  !      do j = 1,NGLLZ
  !        do i = 1,NGLLX
  !          iglob = ibool(i,j,ispec)
  !          if (.not. iglob_is_forced(iglob)) then
  !            ! assembles the contributions
  !            tempx1l = 0._CUSTOM_REAL
  !            tempx2l = 0._CUSTOM_REAL
  !            tempz1l = 0._CUSTOM_REAL
  !            tempz2l = 0._CUSTOM_REAL
  !            do k = 1,NGLLX
  !              tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
  !              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
  !              tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
  !              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
  !            enddo
  !            ! sums contributions from each element to the global values
  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
  !            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)

  !            accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
  !          endif
  !        enddo
  !      enddo
  !    endif
  !  else ! 12.2 lucas
      !lucas7, CTD-SEM
      ! default (not AXISYM case)
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (.not. iglob_is_forced(iglob)) then
            ! assembles the contributions
            tempx1l_m1 = 0._CUSTOM_REAL
            tempx2l_m1 = 0._CUSTOM_REAL
            tempz1l_m1 = 0._CUSTOM_REAL
            tempz2l_m1 = 0._CUSTOM_REAL
            do k = 1,NGLLX
              tempx1l_m1 = tempx1l_m1 + tempx1_m1(k,j) * hprimewgll_xx(k,i)
              tempx2l_m1 = tempx2l_m1 + tempx2_m1(i,k) * hprimewgll_zz(k,j)
              tempz1l_m1 = tempz1l_m1 + tempz1_m1(k,j) * hprimewgll_xx(k,i)
              tempz2l_m1 = tempz2l_m1 + tempz2_m1(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            accel_elastic_m1(1,iglob) = accel_elastic_m1(1,iglob) - (wzgll(j) * tempx1l_m1 + wxgll(i) * tempx2l_m1)
            accel_elastic_m1(2,iglob) = accel_elastic_m1(2,iglob) - (wzgll(j) * tempz1l_m1 + wxgll(i) * tempz2l_m1)
          endif
        enddo
      enddo
  !  endif ! AXISYM

    !13.lucas
    !lucas8, CTD-SEM
    ! adds PML_BOUNDARY_CONDITIONS contribution
 !   if (PML_BOUNDARY_CONDITIONS) then
 !     if (ispec_is_PML(ispec)) then
 !       do j = 1,NGLLZ
 !         do i = 1,NGLLX
 !           iglob = ibool(i,j,ispec)
 !           if (.not. iglob_is_forced(iglob)) then
 !             accel_elastic(1,iglob) = accel_elastic(1,iglob) - accel_elastic_PML(1,i,j)
 !             accel_elastic(2,iglob) = accel_elastic(2,iglob) - accel_elastic_PML(2,i,j)
 !           endif
 !         enddo
 !       enddo
 !     endif
 !   endif

  enddo ! end of loop over all spectral elements

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_4comp_singleA(x1,x2,z1,z2,A,B,C)

! matrix x matrix multiplication, merging 4 loops for x1,x2 = A^t B^t and z1,z2 = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          duz_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!          duz_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprime_xx(i,k)
!            duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + dummy_loc(1,i,k)*hprime_zz(j,k)
!            duz_dgamma(i,j) = duz_dgamma(i,j) + dummy_loc(2,i,k)*hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NDIM,NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x1,x2,z1,z2

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: A
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x1(i,j) = A(1,1,j) * B(i,1) + A(1,2,j) * B(i,2) + A(1,3,j) * B(i,3) + A(1,4,j) * B(i,4) + A(1,5,j) * B(i,5)
        x2(i,j) = A(2,1,j) * B(i,1) + A(2,2,j) * B(i,2) + A(2,3,j) * B(i,3) + A(2,4,j) * B(i,4) + A(2,5,j) * B(i,5)

        z1(i,j) = A(1,i,1) * C(j,1) + A(1,i,2) * C(j,2) + A(1,i,3) * C(j,3) + A(1,i,4) * C(j,4) + A(1,i,5) * C(j,5)
        z2(i,j) = A(2,i,1) * C(j,1) + A(2,i,2) * C(j,2) + A(2,i,3) * C(j,3) + A(2,i,4) * C(j,4) + A(2,i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x1(i,j) = 0._CUSTOM_REAL
        x2(i,j) = 0._CUSTOM_REAL
        z1(i,j) = 0._CUSTOM_REAL
        z2(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x1(i,j) = x1(i,j) + A(1,k,j) * B(i,k)
          x2(i,j) = x2(i,j) + A(2,k,j) * B(i,k)

          z1(i,j) = z1(i,j) + A(1,i,k) * C(j,k)
          z2(i,j) = z2(i,j) + A(2,i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_4comp_singleA

  end subroutine compute_forces_viscoelastic_m1 ! lucas, CTD-SEM


