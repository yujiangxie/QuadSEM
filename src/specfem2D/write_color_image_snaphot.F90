!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently maNZ_IMAGE_color more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

  subroutine write_color_image_snaphot(plot_b_wavefield_only) ! lucas, plot_b_wavefield_only = .false.

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,NGLLX,NGLLZ,REMOVE_PMLS_FROM_JPEG_IMAGES

  use specfem_par, only: myrank,nspec,it,NPROC, &
                        assign_external_model,ibool,kmato,density,rhoext,P_SV, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        displ_elastic,veloc_elastic,accel_elastic, &
                        displ_elastic_m2,displ_elastic_m1,& !lucas, CTD-SEM
                        displs_poroelastic,velocs_poroelastic,accels_poroelastic, &
                        b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                        b_displ_elastic,b_veloc_elastic,b_accel_elastic, &
                        b_displ_elastic_m2, & !lucas, CTD-SEM
                        b_displs_poroelastic,b_velocs_poroelastic,b_accels_poroelastic,SIMULATION_TYPE, &
                        UNDO_ATTENUATION_AND_OR_PML,NO_BACKWARD_RECONSTRUCTION, &
                        CTD_SEM, Full_Hessian_by_Wavefield_Stored,compute_appro_Hessian !lucas, CTD-SEM
                        !no_backward_displ_buffer_lucas_tmp

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  use specfem_par_movie, only: vector_field_display,image_color_data,iglob_image_color, & 
    imagetype_JPEG,nb_pixel_loc, &
    NX_IMAGE_color,NZ_IMAGE_color, &
    num_pixel_loc, &
    vector_field_display_m2, image_color_data_m2,& !lucas, CTD-SEM
    vector_field_display_m1, image_color_data_m1 !lucas, CTD-SEM

#ifdef USE_MPI
  use specfem_par_movie, only: data_pixel_recv,data_pixel_send,nb_pixel_per_proc,num_pixel_recv, &
                               data_pixel_recv_m2, data_pixel_send_m2,& !lucas, CTD-SEM
                               data_pixel_recv_m1, data_pixel_send_m1 !lucas, CTD-SEM
#endif

  implicit none

  !parameter useful for UNDO_ATTENUATION
  logical :: plot_b_wavefield_only

  !local variables
  integer :: i,j,k,ispec,iglob,iproc,i_field,n_fields
  double precision :: rhol

  if (SIMULATION_TYPE == 1 .or. UNDO_ATTENUATION_AND_OR_PML .or. NO_BACKWARD_RECONSTRUCTION) then
    n_fields = 1
  else
    n_fields = 2
  endif

  !loop over each wavefield
  do i_field = 1, n_fields

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color,' for time step ',it
    call flush_IMAIN()
  endif
  
  if (imagetype_JPEG >= 1 .and. imagetype_JPEG <= 3) then !lucas1, CTD-SEM
    if (i_field == 1 .and. .not. plot_b_wavefield_only) then
      if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
        write(IMAIN,*) 'drawing scalar image of the forward wavefield displacement...'
      else if (myrank == 0) then
        write(IMAIN,*) 'drawing scalar image of the adjoint wavefield displacement...'
      endif
      !lucas, from field_elastic to vector_field_element to vector_field_display
      call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic) ! lucas, old version for u*
      !call compute_vector_whole_medium(potential_acoustic,no_backward_displ_buffer_lucas_tmp,displs_poroelastic) ! lucas, changed here for ploting du*
      !call compute_vector_whole_medium(potential_acoustic,b_displ_elastic,displs_poroelastic) ! lucas, changed here for ploting du
      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas, CTD-SEM
      call compute_vector_whole_medium_m2(potential_acoustic,displ_elastic_m2,displs_poroelastic) ! lucas, only implemented for elastic 
      !lucas, from field_elastic_m2 to vector_field_element_m2 to vector_field_display_m2 (see specfem2D_par, prepare_timerun_postscripts())
      !lucas, the forward or adjoint wavefields to be displayed are input here
        if(compute_appro_Hessian) then !for approximate Hessian, 
           call compute_vector_whole_medium_m1(potential_acoustic,displ_elastic_m1,displs_poroelastic) ! lucas, only implemented for elastic 
        endif
      endif
    else
       if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed (or read in) forward wavefield displacement...'
       call compute_vector_whole_medium(b_potential_acoustic,b_displ_elastic,b_displs_poroelastic)
       if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas, CTD-SEM-------
       call compute_vector_whole_medium_m2(b_potential_acoustic,b_displ_elastic_m2,b_displs_poroelastic) ! lucas, only implemented for elastic 
       !lucas, from b_field_elastic_m2 to vector_field_element_m2 to vector_field_display_m2
       !lucas, the backward wavefields to be displayed are input here
       endif!--------     
    endif

  else if (imagetype_JPEG >= 4 .and. imagetype_JPEG <= 6) then

    if (i_field == 1 .and. .not. plot_b_wavefield_only) then
      if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
        write(IMAIN,*) 'drawing scalar image of the forward wavefield velocity...'
      else if (myrank == 0) then
        write(IMAIN,*) 'drawing scalar image of the adjoint wavefield velocity...'
      endif
      call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic)
    else
       if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed forward wavefield velocity...'
      call compute_vector_whole_medium(b_potential_dot_acoustic,b_veloc_elastic,b_velocs_poroelastic)
    endif

  else if (imagetype_JPEG >= 7 .and. imagetype_JPEG <= 9) then

    if (i_field == 1 .and. .not. plot_b_wavefield_only) then
      if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
        write(IMAIN,*) 'drawing scalar image of the forward wavefield acceleration...'
      else if (myrank == 0) then
        write(IMAIN,*) 'drawing scalar image of the adjoint wavefield acceleration...'
      endif
      call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic)
    else
       if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed forward wavefield acceleration...'
      call compute_vector_whole_medium(b_potential_dot_dot_acoustic,b_accel_elastic,b_accels_poroelastic)
    endif

  else if (imagetype_JPEG >= 11 .and. imagetype_JPEG <= 13) then
    ! allocation for normalized representation in JPEG image
    ! for an atmosphere model
    if (myrank == 0) write(IMAIN,*) 'drawing scalar image of part of normalized displacement vector...'
    call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic)

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
          else
            rhol = density(1,kmato(ispec))
          endif
          iglob = ibool(i,j,ispec)
          vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
          vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
        enddo
      enddo
    enddo

  else if (imagetype_JPEG >= 14 .and. imagetype_JPEG <= 16) then
    ! allocation for normalized representation in JPEG image
    ! for an atmosphere model
    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic)

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
          else
            rhol = density(1,kmato(ispec))
          endif
          iglob = ibool(i,j,ispec)
          vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
          vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
        enddo
      enddo
    enddo

  else if (imagetype_JPEG == 10 .and. P_SV) then

    if (i_field == 1 .and. .not. plot_b_wavefield_only) then
      if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
        write(IMAIN,*) 'drawing scalar image of the forward wavefield pressure...'
      else
        if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the adjoint wavefield pressure...'
      endif
    call compute_pressure_whole_medium(1)
    else
      if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed forward wavefield pressure...'
      call compute_pressure_whole_medium(2)
    endif

  else if (imagetype_JPEG == 10 .and. .not. P_SV) then
    call exit_MPI(myrank,'cannot draw pressure field for SH (membrane) waves')

  else
    call exit_MPI(myrank,'wrong type for JPEG snapshots')
  endif

!! DK DK quick hack to remove the PMLs from JPEG images if needed: set the vector field to zero there
  if (PML_BOUNDARY_CONDITIONS .and. REMOVE_PMLS_FROM_JPEG_IMAGES) then  ! lucas: REMOVE_PMLS_FROM_JPEG_IMAGES=true
    do ispec = 1,nspec
      if (ispec_is_PML(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            vector_field_display(1,iglob) = 0.d0
            vector_field_display(2,iglob) = 0.d0
          enddo
        enddo
      endif
    enddo
  endif
!! DK DK quick hack to remove the PMLs from JPEG images if needed

  !lucas2, CTD-SEM
  image_color_data(:,:) = 0.d0
  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
   image_color_data_m2(:,:) = 0.d0
   if(compute_appro_Hessian) then
     image_color_data_m1(:,:) = 0.d0
   endif
  endif

  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

    ! avoid edge effects
    if (i < 1 ) i = 1
    if (j < 1 ) j = 1

    if (i > NX_IMAGE_color ) i = NX_IMAGE_color
    if (j > NZ_IMAGE_color ) j = NZ_IMAGE_color

    if (P_SV) then !lucas3, CTD-SEM
      ! P-SV waves, plot a component of vector, its norm, or else pressure
      if (iglob_image_color(i,j) /= -1) then
        if (imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. &
            imagetype_JPEG == 11 .or. imagetype_JPEG == 14) then
          ! draw the X component of the vector
          image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))

          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
           image_color_data_m2(i,j) = vector_field_display_m2(1,iglob_image_color(i,j))
           if(compute_appro_Hessian) then
             image_color_data_m1(i,j) = vector_field_display_m1(1,iglob_image_color(i,j))
           endif
          endif

        else if (imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. &
                imagetype_JPEG == 12 .or. imagetype_JPEG == 15) then
          ! draw the Z component of the vector
          image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))

          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
           image_color_data_m2(i,j) = vector_field_display_m2(2,iglob_image_color(i,j))
           if(compute_appro_Hessian) then
              image_color_data_m1(i,j) = vector_field_display_m1(2,iglob_image_color(i,j))
           endif
          endif

        else if (imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. &
                imagetype_JPEG == 13 .or. imagetype_JPEG == 16) then
          ! draw the norm of the vector
          image_color_data(i,j) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2  &
                                     + vector_field_display(2,iglob_image_color(i,j))**2)

          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
           image_color_data_m2(i,j) = sqrt(vector_field_display_m2(1,iglob_image_color(i,j))**2  &
                                     + vector_field_display_m2(2,iglob_image_color(i,j))**2)
           if(compute_appro_Hessian) then
              image_color_data_m1(i,j) = sqrt(vector_field_display_m1(1,iglob_image_color(i,j))**2  &
                                     + vector_field_display_m1(2,iglob_image_color(i,j))**2)
           endif
          endif

        else if (imagetype_JPEG == 10) then
          ! by convention we have stored pressure in the 2. component of the array
          image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))
          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
           image_color_data_m2(i,j) = vector_field_display_m2(2,iglob_image_color(i,j))
           if(compute_appro_Hessian) then
              image_color_data_m1(i,j) = vector_field_display_m1(2,iglob_image_color(i,j))
           endif
          endif

        else
          call exit_MPI(myrank,'wrong type for JPEG snapshots')
        endif
      endif

    else
      ! SH (membrane) waves, plot y-component
      if (iglob_image_color(i,j) /= -1) image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))
      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then 
       if (iglob_image_color(i,j) /= -1) image_color_data_m2(i,j) = vector_field_display_m2(1,iglob_image_color(i,j))
       if(compute_appro_Hessian) then
          if (iglob_image_color(i,j) /= -1) image_color_data_m1(i,j) = vector_field_display_m1(1,iglob_image_color(i,j))
       endif
      endif 

    endif
  enddo

  ! assembling array image_color_data on process zero for color output
#ifdef USE_MPI
  if (NPROC > 1) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        call recv_dp(data_pixel_recv(1), nb_pixel_per_proc(iproc), iproc, 43)

        if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas4, CTD-SEM
         call recv_dp(data_pixel_recv_m2(1), nb_pixel_per_proc(iproc), iproc, 430)
         if(compute_appro_Hessian) then
            call recv_dp(data_pixel_recv_m1(1), nb_pixel_per_proc(iproc), iproc, 430+100)
         endif
        endif

        do k = 1, nb_pixel_per_proc(iproc)
          j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
          i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

          ! avoid edge effects
          if (i < 1) i = 1
          if (j < 1) j = 1

          if (i > NX_IMAGE_color) i = NX_IMAGE_color
          if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

          image_color_data(i,j) = data_pixel_recv(k)
          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas5, CTD-SEM
           image_color_data_m2(i,j) = data_pixel_recv_m2(k) 
           if(compute_appro_Hessian) then
             image_color_data_m1(i,j) = data_pixel_recv_m1(k) 
           endif
          endif

        enddo
      enddo
    else !lucas, myrank > 0
      do k = 1, nb_pixel_loc
        j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
        i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

        ! avoid edge effects
        if (i < 1) i = 1
        if (j < 1) j = 1

        if (i > NX_IMAGE_color) i = NX_IMAGE_color
        if (j > NZ_IMAGE_color) j = NZ_IMAGE_color
        !lucas6, CTD-SEM
        if (P_SV) then ! P-SV waves, plot a component of vector, its norm, or else pressure

          if (imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. &
              imagetype_JPEG == 11 .or. imagetype_JPEG == 14) then
             data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector
             if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
              data_pixel_send_m2(k) = vector_field_display_m2(1,iglob_image_color(i,j))
              if(compute_appro_Hessian) then
                data_pixel_send_m1(k) = vector_field_display_m1(1,iglob_image_color(i,j))
              endif
             endif

          else if (imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. &
                  imagetype_JPEG == 12 .or. imagetype_JPEG == 15) then
             data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))  ! draw the Z component of the vector
             if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
              data_pixel_send_m2(k) = vector_field_display_m2(2,iglob_image_color(i,j))
              if(compute_appro_Hessian)then
                data_pixel_send_m1(k) = vector_field_display_m1(2,iglob_image_color(i,j))
              endif
             endif
          else if (imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. &
                  imagetype_JPEG == 13 .or. imagetype_JPEG == 16) then
            data_pixel_send(k) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                      vector_field_display(2,iglob_image_color(i,j))**2)  ! draw the norm of the vector
            if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
             data_pixel_send_m2(k) = sqrt(vector_field_display_m2(1,iglob_image_color(i,j))**2 + &
                                      vector_field_display_m2(2,iglob_image_color(i,j))**2)  ! draw the norm of the vector
             if(compute_appro_Hessian) then
               data_pixel_send_m1(k) = sqrt(vector_field_display_m1(1,iglob_image_color(i,j))**2 + &
                                      vector_field_display_m1(2,iglob_image_color(i,j))**2)  ! draw the norm of the vector
             endif
            endif

          else if (imagetype_JPEG == 10) then
            ! by convention we have stored pressure in the 2. component of the array
            data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))

          else
            call exit_MPI(myrank,'wrong type for JPEG snapshots')
          endif

        else ! SH (membrane) waves, plot y-component
          if (iglob_image_color(i,j) /= -1) data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))
          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
           if (iglob_image_color(i,j) /= -1) data_pixel_send_m2(k) = vector_field_display_m2(1,iglob_image_color(i,j))
           if(compute_appro_Hessian) then
             if (iglob_image_color(i,j) /= -1) data_pixel_send_m1(k) = vector_field_display_m1(1,iglob_image_color(i,j))
           endif
          endif
        endif
      enddo !lucas, 

      call send_dp(data_pixel_send(1), nb_pixel_loc, 0, 43)
      !lucas7, CTD-SEM
      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then
       call send_dp(data_pixel_send_m2(1), nb_pixel_loc, 0, 430)
       if(compute_appro_Hessian) then
        call send_dp(data_pixel_send_m1(1), nb_pixel_loc, 0, 430+100)
       endif
      endif

    endif ! end myrank
  endif
  call synchronize_all()
#else
  ! dummy to avoid compiler warning
  iproc = NPROC
#endif

  call synchronize_all()

  ! creates image
  if (myrank == 0) then
    call create_color_image(i_field,plot_b_wavefield_only) 

    ! user output
    write(IMAIN,*) 'Color image created'
    call flush_IMAIN()
  endif
  call synchronize_all()

enddo !loop over wavefields (forward and adjoint if necessary)

  end subroutine write_color_image_snaphot

