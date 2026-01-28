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

  subroutine finalize_simulation()

#ifdef USE_MPI
  use mpi
#endif
  use constants, only: TWO,FOUR_THIRDS,TWO_THIRDS,APPROXIMATE_HESS_KL,IMAIN,IOUT_ENERGY,ISTANDARD_OUTPUT,IN_DATA_FILES,OUTPUT_FILES
  use specfem_par
  use specfem_par_noise
  use specfem_par_gpu
  use specfem_par_movie, only: simulation_title,output_wavefield_dumps,mask_ibool

  implicit none

  integer :: i,ispec,j,iglob
  integer :: ier
  real(kind=4),dimension(:,:,:),allocatable :: rho_save, vp_save, vs_save, kappa_save, x_save, z_save, Qkappa_save,Qmu_save
  real(kind=4),dimension(:,:,:),allocatable :: rho_save_m2, vp_save_m2,vs_save_m2, Qkappa_save_m2,Qmu_save_m2 ! CTD-SEM
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic
  character(len=MAX_STRING_LEN) :: inputname,outputname,outputname2
  ! 1. lucas
  ! writes out kernel files
  if (SIMULATION_TYPE == 3) then
    call save_adjoint_kernels() !lucas, CTD-SEM, alread done for elastic case.
  endif

  ! saves model files
  if (trim(SAVE_MODEL) /= 'default' .and. trim(SAVE_MODEL) /= '.false.') then

    ! allocates temporary arrays for file storage
    allocate(rho_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 01')

    allocate(vp_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 02')

    allocate(vs_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 03')

    allocate(kappa_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 04')

    allocate(x_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 05')

    allocate(z_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 06')

    if (ATTENUATION_VISCOACOUSTIC) then
      allocate(Qkappa_save(NGLLX,NGLLZ,nspec),stat=ier)
      if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 07')
    endif

    if (ATTENUATION_VISCOELASTIC) then
      if (ATTENUATION_VISCOACOUSTIC) call exit_MPI(myrank, &
                    'Not possible yet to save model with both acoustic and elastic attenuation')
      allocate(Qkappa_save(NGLLX,NGLLZ,nspec),Qmu_save(NGLLX,NGLLZ,nspec),stat=ier)
      if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 08')
    endif

    do ispec= 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          rho_save(i,j,ispec) = density(1,kmato(ispec))
          lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
          mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
          if (AXISYM) then ! CHECK kappa
            kappa_save(i,j,ispec) = lambdal_unrelaxed_elastic + TWO_THIRDS * mul_unrelaxed_elastic
            vp_save(i,j,ispec) = sqrt((kappa_save(i,j,ispec) + FOUR_THIRDS *mul_unrelaxed_elastic)/density(1,kmato(ispec)))
          else
            kappa_save(i,j,ispec) = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
            vp_save(i,j,ispec) = sqrt((kappa_save(i,j,ispec) + mul_unrelaxed_elastic)/density(1,kmato(ispec)))
          endif

          vs_save(i,j,ispec) = sqrt(mul_unrelaxed_elastic/rho_save(i,j,ispec))

          if (assign_external_model) then
            vp_save(i,j,ispec) = vpext(i,j,ispec)
            vs_save(i,j,ispec) = vsext(i,j,ispec)
            rho_save(i,j,ispec) = rhoext(i,j,ispec)
          endif

          iglob = ibool(i,j,ispec)
          x_save(i,j,ispec) = coord(1,iglob)
          z_save(i,j,ispec) = coord(2,iglob)
          if (ATTENUATION_VISCOACOUSTIC) Qkappa_save(i,j,ispec) = QKappa_attenuation(kmato(ispec))

          if (ATTENUATION_VISCOELASTIC) then !
            Qkappa_save(i,j,ispec) = QKappa_attenuation(kmato(ispec))
            Qmu_save(i,j,ispec) = Qmu_attenuation(kmato(ispec))

            if (assign_external_model) then !------only added for the visco-ela
              Qkappa_save(i,j,ispec) = QKappa_attenuationext(i,j,ispec)
              Qmu_save(i,j,ispec) = Qmu_attenuationext(i,j,ispec)
            endif !-----------------------------------------------------

          endif
        enddo
      enddo
    enddo

   !lucas add this for m2-------------------------------------------------------------------------------------------
    if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
      ! allocates temporary arrays for file storage
      allocate(rho_save_m2(NGLLX,NGLLZ,nspec),stat=ier)
      if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 01_m2')

      allocate(vp_save_m2(NGLLX,NGLLZ,nspec),stat=ier)
      if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 02_m2')

      allocate(vs_save_m2(NGLLX,NGLLZ,nspec),stat=ier)
      if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 03_m2')

      if (ATTENUATION_VISCOACOUSTIC) then
        allocate(Qkappa_save_m2(NGLLX,NGLLZ,nspec),stat=ier)
        if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 07_m2')
      endif

      if (ATTENUATION_VISCOELASTIC) then
        if (ATTENUATION_VISCOACOUSTIC) call exit_MPI(myrank, &
                    'Not possible yet to save model with both acoustic and elastic attenuation')
        allocate(Qkappa_save_m2(NGLLX,NGLLZ,nspec),Qmu_save_m2(NGLLX,NGLLZ,nspec),stat=ier)
        if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 08_m2')
      endif

      do ispec= 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            if (assign_external_model) then
              vp_save_m2(i,j,ispec) = vpext_m2(i,j,ispec)
              vs_save_m2(i,j,ispec) = vsext_m2(i,j,ispec)
              rho_save_m2(i,j,ispec) = rhoext_m2(i,j,ispec)
            endif

            if (ATTENUATION_VISCOACOUSTIC) Qkappa_save_m2(i,j,ispec) = QKappa_attenuationext_m2(i,j,ispec)
            if (ATTENUATION_VISCOELASTIC) then
              Qkappa_save_m2(i,j,ispec) = QKappa_attenuationext_m2(i,j,ispec)
              Qmu_save_m2(i,j,ispec) = Qmu_attenuationext_m2(i,j,ispec)
            endif
          enddo
        enddo
      enddo

    endif ! CTD-SEM end 
   !----------------------------------------------------------------------------------------------------------------


    ! SMNSR For compatibility with NUMBER_OF_SIMULTANEOUS_RUNS we have to change the lines trim(IN_DATA_FILES)//'proc'

    ! outputs model files
    if (trim(SAVE_MODEL) == 'legacy') then
      ! legacy format
      write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_model_velocity.dat_input'
      open(unit=1001,file=inputname,status='unknown',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_model_velocity.dat_input')
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            write(1001,'(I10,5e15.5e4)') iglob, x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec), &
                                     vp_save(i,j,ispec),vs_save(i,j,ispec)
          enddo
        enddo
      enddo
      close(1001)

    else if (trim(SAVE_MODEL) == 'ascii') then
      ! ascii format
      write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho_vp_vs.dat'
      open(unit=1001,file= inputname,status='unknown',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_rho_vp_vs.dat')
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            write(1001,'(5e15.5e4)') x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec),vp_save(i,j,ispec),vs_save(i,j,ispec)
          enddo
        enddo
      enddo
      close(1001)

    else if ((trim(SAVE_MODEL) == 'binary') .or. (trim(SAVE_MODEL) == 'gll')) then
      ! binary and GLL format
      write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_rho.bin')
      write(172) rho_save
      close(172)
      write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vp.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vp.bin')
      write(172) vp_save
      close(172)
      write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vs.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vs.bin')
      write(172) vs_save
      close(172)
      write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_x.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_x.bin')
      write(172) x_save
      close(172)
      write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_z.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_z.bin')
      write(172) z_save
      close(172)

      write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_jacobian.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_jacobian.bin')
      write(172) jacobian
      close(172)

      if (ATTENUATION_VISCOACOUSTIC) then
        write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa.bin'
        open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qkappa.bin')
        write(172) Qkappa_save
        close(172)
      endif

      if (ATTENUATION_VISCOELASTIC) then
        write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa.bin'
        open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qkappa.bin')
        write(172) Qkappa_save
        close(172)

        write(outputname,'(a,i6.6,a)')trim(IN_DATA_FILES)//'proc',myrank,'_Qmu.bin'
        open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qmu.bin')
        write(172) Qmu_save
        close(172)
      endif

      ! lucas added this for saving gll models for m2---------------------------------------
      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas, CTD-SEM
        ! binary and GLL format
        write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho_m2.bin'
        open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_rho_m2.bin')
        write(172) rho_save_m2
        close(172)
        write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vp_m2.bin'
        open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vp_m2.bin')
        write(172) vp_save_m2
        close(172)
        write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vs_m2.bin'
        open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vs_m2.bin')
        write(172) vs_save_m2
        close(172)
      !  write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_x.bin'
      !  open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      !  if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_x.bin')
      !  write(172) x_save
      !  close(172)
      !  write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_z.bin'
      !  open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      !  if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_z.bin')
      !  write(172) z_save
      !  close(172)

      !  write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_jacobian.bin'
      !  open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      !  if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_jacobian.bin')
      !  write(172) jacobian
      !  close(172)

        if (ATTENUATION_VISCOACOUSTIC) then
          write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa_m2.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qkappa_m2.bin')
          write(172) Qkappa_save_m2
          close(172)
        endif

        if (ATTENUATION_VISCOELASTIC) then
          write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa_m2.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qkappa_m2.bin')
          write(172) Qkappa_save_m2
          close(172)

          write(outputname,'(a,i6.6,a)')trim(IN_DATA_FILES)//'proc',myrank,'_Qmu_m2.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qmu_m2.bin')
          write(172) Qmu_save_m2
          close(172)
        endif

      endif ! end for CTD-SEM
      !-------------------------------------------------------------------------------------


    else
      call stop_the_code('Save Model not implemented for external and tomo')
    endif !Type of model
  endif !save model


  ! 2. lucas 
  ! For this mode, the forward model has been saved differently
  if (((.not. NO_BACKWARD_RECONSTRUCTION) .or. NO_BACKWARD_RECONSTRUCTION) .and. (.not. ATTENUATION_VISCOELASTIC)) then !lucas, add ATTENUATION_VISCOELASTIC, need to check when PML

    ! stores absorbing boundary contributions into files
    if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. STACEY_ABSORBING_CONDITIONS) then

      if (any_acoustic) then
        !--- left absorbing boundary
        if (nspec_left > 0) then
          ! opens file
          write(outputname,'(a,i6.6,a)') 'absorb_acoustic_left',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_acoustic_left
          close(35)
        endif
        !--- right absorbing boundary
        if (nspec_right > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_acoustic_right',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_acoustic_right
          close(35)
        endif
        !--- bottom absorbing boundary
        if (nspec_bottom > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_acoustic_bottom',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_acoustic_bottom
          close(35)
        endif
        !--- top absorbing boundary
        if (nspec_top > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_acoustic_top',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_acoustic_top
          close(35)
        endif
      endif !any acoustic

      if (any_elastic) then
        !--- left absorbing boundary
        if (nspec_left > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_left',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_elastic_left
          close(35)
        endif
        !--- right absorbing boundary
        if (nspec_right > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_right',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_elastic_right
          close(35)
        endif
        !--- bottom absorbing boundary
        if (nspec_bottom > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_bottom',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_elastic_bottom
          close(35)
        endif
        !--- top absorbing boundary
        if (nspec_top > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_top',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_elastic_top
          close(35)
        endif

      endif !any elastic

     ! lucas, CTD-SEM ------------------- in prepare_absorb.f90 the open for read should be unit=350, not 35, lucas-------------------------
      if(CTD_SEM) then !lucas, save the absorbing boundary fields
      if (any_elastic) then
        !--- left absorbing boundary
        if (nspec_left > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_left_m2_',myrank,'.bin'
          open(unit=350,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(350) b_absorb_elastic_left_m2
          close(350)
        endif
        !--- right absorbing boundary
        if (nspec_right > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_right_m2_',myrank,'.bin'
          open(unit=350,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(350) b_absorb_elastic_right_m2
          close(350)
        endif
        !--- bottom absorbing boundary
        if (nspec_bottom > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_bottom_m2_',myrank,'.bin'
          open(unit=350,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(350) b_absorb_elastic_bottom_m2
          close(350)
        endif
        !--- top absorbing boundary
        if (nspec_top > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_elastic_top_m2_',myrank,'.bin'
          open(unit=350,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(350) b_absorb_elastic_top_m2
          close(350)
        endif
      endif !any elastic
      endif
      ! --------------------------------------------------------------------------------------------------

      if (any_poroelastic) then
        !--- left absorbing boundary
        if (nspec_left > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_poro_s_left',myrank,'.bin'
          write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_left',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_poro_s_left
          write(36) b_absorb_poro_w_left
          close(35)
          close(36)
        endif
        !--- right absorbing boundary
        if (nspec_right > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_poro_s_right',myrank,'.bin'
          write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_right',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_poro_s_right
          write(36) b_absorb_poro_w_right
          close(35)
          close(36)
        endif
        !--- bottom absorbing boundary
        if (nspec_bottom > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_poro_s_bottom',myrank,'.bin'
          write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_bottom',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_poro_s_bottom
          write(36) b_absorb_poro_w_bottom
          close(35)
          close(36)
        endif
        !--- top absorbing boundary
        if (nspec_top > 0) then
          write(outputname,'(a,i6.6,a)') 'absorb_poro_s_top',myrank,'.bin'
          write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_top',myrank,'.bin'
          open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
          ! writes boundary contributions
          write(35) b_absorb_poro_s_top
          write(36) b_absorb_poro_w_top
          close(35)
          close(36)
        endif

      endif !lucas, any_poroelastic
    endif ! lucas, SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. STACEY_ABSORBING_CONDITIONS

    ! PML
    if (anyabs_glob .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. PML_BOUNDARY_CONDITIONS) then
      if (any_elastic .and. nglob_interface > 0) close(71)
      if (any_acoustic .and. nglob_interface > 0) close(72)
    endif

    ! save last frame
    !------------------------------- elastic-------------------------------------------------------------
    if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_elastic) then
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Saving elastic last frame'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
     ! (1)
      write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
      open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_elastic**.bin')

      write(55) displ_elastic
      write(55) veloc_elastic
      write(55) accel_elastic
      close(55)
    endif
    ! lucas, CTD-SEM--------------------------------------------------------------------------------------
     if(CTD_SEM) then ! lucas, save the last-state field
     if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_elastic) then
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Saving elastic last frame for m2 in CTD-SEM'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
     ! (1)
      write(outputname,'(a,i6.6,a)') 'lastframe_elastic_m2_',myrank,'.bin'
      open(unit=550,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_elastic_m2_**.bin for m2 write in CTD-SEM')

      write(550) displ_elastic_m2
      write(550) veloc_elastic_m2
      write(550) accel_elastic_m2
      close(550)
    endif
    endif
    
    !----------------------------------------------------------------------------------------------------

    if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_poroelastic) then
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Saving poroelastic last frame...'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
      open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_poroelastic_s**.bin')

      write(55) displs_poroelastic
      write(55) velocs_poroelastic
      write(55) accels_poroelastic
      close(55)

      write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
      open(unit=56,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_poroelastic_w**.bin')

      write(56) displw_poroelastic
      write(56) velocw_poroelastic
      write(56) accelw_poroelastic
      close(56)
    endif

    if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_acoustic) then
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Saving acoustic last frame...'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
      open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_acoustic**.bin')

      write(55) potential_acoustic
      write(55) potential_dot_acoustic
      write(55) potential_dot_dot_acoustic
      close(55)
    endif

  endif ! if trim(SAVE_MODEL) /= '.false.' .and. (.not.
        ! UNDO_ATTENUATION_AND_OR_PML) .and. (.not. NO_BACKWARD_RECONSTRUCTION)
  ! 3.lucas 
  ! frees memory
  if (GPU_MODE) then
    ! frees temporary arrays
    if (any_elastic) deallocate(tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D)

    deallocate(request_send_recv_scalar_gpu,b_request_send_recv_scalar_gpu)
    deallocate(request_send_recv_vector_gpu,b_request_send_recv_vector_gpu)
    deallocate(buffer_send_scalar_gpu,b_buffer_send_scalar_gpu)
    deallocate(buffer_recv_scalar_gpu,b_buffer_recv_scalar_gpu)
    deallocate(buffer_send_vector_gpu,b_buffer_send_vector_gpu)
    deallocate(buffer_recv_vector_gpu,b_buffer_recv_vector_gpu)

    ! frees memory on GPU
    call prepare_cleanup_device(Mesh_pointer,any_acoustic,any_elastic, &
                                ANISOTROPY, &
                                APPROXIMATE_HESS_KL, &
                                ATTENUATION_VISCOACOUSTIC, &
                                ATTENUATION_VISCOELASTIC, &
                                NO_BACKWARD_RECONSTRUCTION, &
                                no_backward_acoustic_buffer)
  endif

  if (output_wavefield_dumps) deallocate(mask_ibool)

  if (initialfield .and. over_critical_angle) then
    deallocate(v0x_left)
    deallocate(v0z_left)
    deallocate(t0x_left)
    deallocate(t0z_left)

    deallocate(v0x_right)
    deallocate(v0z_right)
    deallocate(t0x_right)
    deallocate(t0z_right)

    deallocate(v0x_bot)
    deallocate(v0z_bot)
    deallocate(t0x_bot)
    deallocate(t0z_bot)
  endif
  ! frees arrays (not complete, but at least a few...)
  deallocate(seismo_current)
  deallocate(sisux,sisuz,siscurl)

  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas, CTD-SEM
  deallocate(sisux_m2,sisuz_m2,siscurl_m2)
  endif

  ! wavefields
  deallocate(displ_elastic,veloc_elastic,accel_elastic)
  deallocate(displ_elastic_old)

  !lucas, CTD-SEM-------------
  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
  deallocate(displ_elastic_m2,veloc_elastic_m2,accel_elastic_m2)
  deallocate(displ_elastic_old_m2)
   if(compute_appro_Hessian) then !----
     deallocate(displ_elastic_m1,veloc_elastic_m1,accel_elastic_m1)
     deallocate(displ_elastic_old_m1)
   endif!----
  endif
  !--------------------
  deallocate(rmass_inverse_elastic)

  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
  deallocate(rmass_inverse_elastic_m2) !lucas, CTD-SEM
  endif

  deallocate(b_displ_elastic,b_veloc_elastic,b_accel_elastic)
  deallocate(b_displ_elastic_old)
  
  !CTD-SEM------
  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
  deallocate(b_displ_elastic_m2,b_veloc_elastic_m2,b_accel_elastic_m2)
  deallocate(b_displ_elastic_old_m2)
  endif
  !---------------------

  deallocate(potential_acoustic,potential_acoustic_old)
  deallocate(potential_dot_acoustic,potential_dot_dot_acoustic)
  deallocate(rmass_inverse_acoustic)
  deallocate(b_potential_acoustic,b_potential_acoustic_old)
  if (.not. NO_BACKWARD_RECONSTRUCTION) then
    deallocate(b_potential_dot_acoustic,b_potential_dot_dot_acoustic)
  endif
  deallocate(b_displ_ac,b_accel_ac,accel_ac)
  ! noise
  if (allocated(time_function_noise)) deallocate(time_function_noise)
  if (allocated(source_array_noise)) deallocate(source_array_noise)
  if (allocated(mask_noise)) deallocate(mask_noise)
  if (allocated(surface_movie_y_or_z_noise)) deallocate(surface_movie_y_or_z_noise)
  ! material
  deallocate(kappastore,mustore,rhostore,rho_vp,rho_vs)
  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
  deallocate(kappastore_m2,mustore_m2,rhostore_m2,rho_vp_m2,rho_vs_m2) !lucas, CTD-SEM
  endif
  ! attenuation
  deallocate(e1,e11,e13)
  deallocate(dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old)
  deallocate(A_newmark_nu1,B_newmark_nu1,A_newmark_nu2,B_newmark_nu2)
  deallocate(e1_acous_sf,sum_forces_old)
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      deallocate(b_e1_acous_sf,b_sum_forces_old)
    endif
    if (any_elastic) then
      deallocate(b_e1,b_e11,b_e13,b_dux_dxl_old,b_duz_dzl_old,b_dux_dzl_plus_duz_dxl_old)
    endif
  endif

  ! close energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)

  ! print exit banner
  if (myrank == 0) call datim(simulation_title)

  ! close output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  end subroutine finalize_simulation
