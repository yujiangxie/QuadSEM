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

  subroutine save_forward_arrays_undoatt() !lucas, done also for m2

  use constants, only: IOUT_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,NGLLX,NGLLZ

  use specfem_par, only: myrank,iteration_on_subset, &
    any_acoustic,any_elastic,ATTENUATION_VISCOACOUSTIC,ATTENUATION_VISCOELASTIC, &
    potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
    displ_elastic,veloc_elastic,accel_elastic, &
    e1,e11,e13,dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old, &
    displ_elastic_m2,veloc_elastic_m2,accel_elastic_m2, & !lucas, CTD-SEM
    e1_m2,e11_m2,e13_m2,dux_dxl_old_m2,duz_dzl_old_m2,dux_dzl_plus_duz_dxl_old_m2, CTD_SEM, & !lucas, CTD-SEM
    e1_acous_sf,sum_forces_old,GPU_MODE,nspec_ATT_ac,nglob

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  ! saves frame of the forward simulation

  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
  open(unit=IOUT_UNDO_ATT  ,file=trim(OUTPUT_FILES)//outputname, &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for writing')
 
  if(CTD_SEM) then ! lucas add output for m2 -------------
  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_m2_save_frame_at',iteration_on_subset_tmp,'.bin'
  open(unit=IOUT_UNDO_ATT+1  ,file=trim(OUTPUT_FILES)//outputname, &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_m2_save_frame_at** for writing')
  endif
  !---------------------

  if (any_acoustic) then
    if (GPU_MODE) call transfer_fields_ac_from_device(nglob,potential_acoustic,potential_dot_acoustic, &
                                                      potential_dot_dot_acoustic,Mesh_pointer)
    write(IOUT_UNDO_ATT) potential_dot_dot_acoustic
    write(IOUT_UNDO_ATT) potential_dot_acoustic
    write(IOUT_UNDO_ATT) potential_acoustic

    if (ATTENUATION_VISCOACOUSTIC) then
      if (GPU_MODE) call transfer_viscoacoustic_var_from_device(NGLLX*NGLLZ*nspec_ATT_ac, &
                                                                e1_acous_sf,sum_forces_old,Mesh_pointer)
      write(IOUT_UNDO_ATT) e1_acous_sf
      write(IOUT_UNDO_ATT) sum_forces_old

    endif

  endif

  if (any_elastic) then
    write(IOUT_UNDO_ATT) accel_elastic
    write(IOUT_UNDO_ATT) veloc_elastic
    write(IOUT_UNDO_ATT) displ_elastic

    if (ATTENUATION_VISCOELASTIC) then
      write(IOUT_UNDO_ATT) e1
      write(IOUT_UNDO_ATT) e11
      write(IOUT_UNDO_ATT) e13
      write(IOUT_UNDO_ATT) dux_dxl_old
      write(IOUT_UNDO_ATT) duz_dzl_old
      write(IOUT_UNDO_ATT) dux_dzl_plus_duz_dxl_old
    endif
  endif
  close(IOUT_UNDO_ATT)

  if(CTD_SEM) then ! lucas add for --------------m2
  if (any_elastic) then
    write(IOUT_UNDO_ATT+1) accel_elastic_m2
    write(IOUT_UNDO_ATT+1) veloc_elastic_m2
    write(IOUT_UNDO_ATT+1) displ_elastic_m2

    if (ATTENUATION_VISCOELASTIC) then
      write(IOUT_UNDO_ATT+1) e1_m2
      write(IOUT_UNDO_ATT+1) e11_m2
      write(IOUT_UNDO_ATT+1) e13_m2
      write(IOUT_UNDO_ATT+1) dux_dxl_old_m2
      write(IOUT_UNDO_ATT+1) duz_dzl_old_m2
      write(IOUT_UNDO_ATT+1) dux_dzl_plus_duz_dxl_old_m2
    endif
  endif
  close(IOUT_UNDO_ATT+1)
  endif !--------------------------------------- m2 

  end subroutine save_forward_arrays_undoatt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine save_forward_arrays_no_backward()

  use constants, only: IOUT_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,APPROXIMATE_HESS_KL,NDIM

  use specfem_par, only: myrank, islice_selected_rec,ispec_selected_rec,ibool, it, & !it, NSTEP, &  
    any_acoustic,any_elastic,potential_acoustic,displ_elastic,GPU_MODE, accel_elastic, accel_elastic_m2, &  ! have removed  b_displ_elastic,b_displ_elastic_m2
    displ_elastic_m2, Full_Hessian_by_Wavefield_Stored, & !lucas, CTD-SEM 
    no_backward_acoustic_buffer,no_backward_displ_buffer, & !no_backward_accel_buffer, &  lucas
    no_backward_iframe,no_backward_Nframes,nglob,SIMULATION_TYPE, nspec, coord ! lucas added SIMULATION_TYPE,nspec, coord

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! EB EB June 2018 : in this routine, in order to overlap both GPU = => RAM and RAM = => disk transfers, we transfer
  ! a wavefield in two iterations.
  ! At the first iteration, we transfer the wavefield from the GPU to the disk.
  ! At the second iteration, we write this wavefield on the disk from the RAM.
  ! In the text above, an iteration means NSTEP_BETWEEN_COMPUTE_KERNELS iterations of the timeloop.
  ! The buffer no_backward_acoustic_buffer is declared in only one dimension in
  ! order to allow the CUDA API to set it in pinned memory (HostRegister).
  ! To perform the async I/O, stream accesses are used for files, numerical
  ! experiences showed that it is the fastest way.

  ! local parameters
  integer :: ier,buffer_num_GPU_transfer,buffer_num_async_IO
  integer(KIND=8) :: offset
  character(len=MAX_STRING_LEN) :: outputname
  integer :: iglob_source ! lucas
  integer :: i, j, ispec, iglob !Lucas
  double precision :: xx, zz !Lucas
  logical :: write_field_at_frames, for_m1
  write_field_at_frames=.false. !lucas, need to set in Par_file
  for_m1=.false.
  !logical :: file_exists !lucas

  ! safety check
  if (GPU_MODE .and. any_elastic) call stop_the_code('No backward simulation is not available for elastic on GPU')

  ! increments counter of wavefield frames transfered
  no_backward_iframe = no_backward_iframe + 1   ! lucas: this routine in manage_no_backward_reconstruction_io() of the time loop , no_backward_iframe=0 in prapare_time_run

  ! opens file to save at the first use of the routine
  !if (no_backward_iframe == 1) then
  !  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_No_backward_reconstruction_database.bin'
  ! ! lucas-----------------------------------------------------------------------------------------------------------------
  ! !INQUIRE(FILE=trim(OUTPUT_FILES)//outputname, EXIST=file_exists)
  !    open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
  !       status='unknown',asynchronous='yes',form='unformatted',action='write',access='stream',iostat=ier) 
  !    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_No_backward_reconstruction_database.bin for writing')
  !else if (no_backward_iframe > 2) then
  !  ! waits for previous asynchronous I/O
  !  wait(IOUT_UNDO_ATT)
  !endif
  ! lucas ------------------------------------------------------------------------------------------------------------------
  if(Full_Hessian_by_Wavefield_Stored) then
   if(SIMULATION_TYPE==1) then ! lucas: m1 for u(m) and m2 for u(m+dm)
      
      write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m1_iframe_at',no_backward_iframe,'.bin'
      open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',form='unformatted',action='write',iostat=ier) 
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_m1_iframe_at** for writing')
          if(write_field_at_frames .and. for_m1) then
          !---------
          if(no_backward_iframe==000101 .or. no_backward_iframe==000301 .or. no_backward_iframe==000501  &
          .or. no_backward_iframe==000701 .or. no_backward_iframe==000901) then ! write one time frame for ploting, 1000, 3000, 5000, 7000, 9000 steps
          write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m1_iframe_at',no_backward_iframe,'.dat'
          open(unit=no_backward_iframe,file=trim(OUTPUT_FILES)//outputname, &
          status='unknown',iostat=ier) 
          if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_m1_iframe_at** for .dat writing')
          endif
          !--------
          endif
           
     write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m2_iframe_at',no_backward_iframe,'.bin'
      open(unit=IOUT_UNDO_ATT+100,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',form='unformatted',action='write',iostat=ier) 
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_m2_iframe_at** for writing')
          if(write_field_at_frames .and. (.not. for_m1)) then
          !---------
          if(no_backward_iframe==000101 .or. no_backward_iframe==000301 .or. no_backward_iframe==000501  &
          .or. no_backward_iframe==000701 .or. no_backward_iframe==000901) then ! write one time frame for ploting,
          write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m2_iframe_at',no_backward_iframe,'.dat'
          open(unit=no_backward_iframe,file=trim(OUTPUT_FILES)//outputname, &
          status='unknown',iostat=ier) 
          if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_m2_iframe_at** for .dat writing')
          endif
          !--------
          endif

      ! write for accelaration for computing density kernels in Ha and Hb
      write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_accel_m1_iframe_at',no_backward_iframe,'.bin'
      open(unit=IOUT_UNDO_ATT+200,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',form='unformatted',action='write',iostat=ier) 
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_accel_m1_iframe_at** for writing')
          
      write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_accel_m2_iframe_at',no_backward_iframe,'.bin'
      open(unit=IOUT_UNDO_ATT+300,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',form='unformatted',action='write',iostat=ier) 
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_accel_m2_iframe_at** for writing')


   else if(SIMULATION_TYPE==3) then !lucas, save for adjoint field
     
    !    write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m1_iframe_at',no_backward_iframe,'.bin'
    !    open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
    !    status='unknown',form='unformatted',action='write',iostat=ier) 
    !    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m1_iframe_at** for writing')
          if(write_field_at_frames .and. for_m1) then
          !---------
          if(no_backward_iframe==000101 .or. no_backward_iframe==000301 .or. no_backward_iframe==000501  &
          .or. no_backward_iframe==000701 .or. no_backward_iframe==000901) then ! write one time frame for ploting,
          write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m1_iframe_at',no_backward_iframe,'.dat'
          open(unit=no_backward_iframe,file=trim(OUTPUT_FILES)//outputname, &
          status='unknown',iostat=ier) 
          if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m1_iframe_at** for .dat writing')
          endif
          !--------
          endif

        write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m2_iframe_at',no_backward_iframe,'.bin'
        open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
        status='unknown',form='unformatted',action='write',iostat=ier) 
        if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m2_iframe_at** for writing')
             if(write_field_at_frames .and. (.not. for_m1)) then
             !---------
             if(no_backward_iframe==000101 .or. no_backward_iframe==000301 .or. no_backward_iframe==000501  &
             .or. no_backward_iframe==000701 .or. no_backward_iframe==000901) then ! write one time frame for ploting,
             write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m2_iframe_at',no_backward_iframe,'.dat'
             open(unit=no_backward_iframe,file=trim(OUTPUT_FILES)//outputname, &
             status='unknown',iostat=ier) 
             if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m2_iframe_at** for .dat writing')
             endif
             !--------
             endif
     
   else
      call exit_MPI(myrank,'Error SIMULATION_TYPE for computing the full Hessian kernels')
   endif
 else if (.not. Full_Hessian_by_Wavefield_Stored) then ! the original version
    if(SIMULATION_TYPE==1) then 
      write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m1_iframe_at',no_backward_iframe,'.bin' ! or m2
      open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
           status='unknown',form='unformatted',action='write',iostat=ier) 
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_iframe_at** for writing')

      write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_accel_m1_iframe_at',no_backward_iframe,'.bin'  ! or m2
      open(unit=IOUT_UNDO_ATT+200,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',form='unformatted',action='write',iostat=ier) 
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_accel_m1_iframe_at** for writing')
          
         
    else if(SIMULATION_TYPE==3) then !lucas, save for adjoint field
      write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m1_iframe_at',no_backward_iframe,'.bin'
      open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
           status='unknown',form='unformatted',action='write',iostat=ier) 
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_iframe_at** for writing')
    endif

 else
   call exit_MPI(myrank,'Error SIMULATION_TYPE')
 
 endif ! end Full_Hessian_by_Wavefield_Stored

  buffer_num_GPU_transfer = mod(no_backward_iframe+2,3)
  buffer_num_async_IO = mod(no_backward_iframe,3)

  ! for the two first times, we only launch GPU = => RAM transfers
  !if (no_backward_iframe < 3) then

  !  if (GPU_MODE) then
  !    call transfer_async_pot_ac_from_device(no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1),Mesh_pointer)
  !  else
      no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1:nglob*(buffer_num_GPU_transfer+1)) = potential_acoustic
  !  endif

  !else if (no_backward_iframe <= no_backward_Nframes) then

    if (any_acoustic) then

      ! the offset is calculated in two steps in order to avoid integer overflow
      offset = no_backward_iframe - 3
      offset = offset * nglob * 4 + 1
      write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset) &
                          no_backward_acoustic_buffer(nglob*buffer_num_async_IO+1:nglob*(buffer_num_async_IO+1))

      if (GPU_MODE) then
        call transfer_async_pot_ac_from_device(no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1),Mesh_pointer)
      else
        no_backward_acoustic_buffer(nglob*mod(no_backward_iframe+1,3)+1:nglob*(mod(no_backward_iframe+1,3)+1)) = potential_acoustic
      endif

      ! for the last transfer, we need to add a statement to wait for the last frame
      if (no_backward_iframe == no_backward_Nframes) then
        ! call to finalize disk writing
        wait(IOUT_UNDO_ATT)

        ! call to finalize GPU transfer, which also initiate a (dummy) transfer
        if (GPU_MODE) &
          call transfer_async_pot_ac_from_device(no_backward_acoustic_buffer(nglob*buffer_num_async_IO+1),Mesh_pointer)

        write(IOUT_UNDO_ATT,pos=offset+4*nglob) &
                            no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1:nglob*(buffer_num_GPU_transfer+1))
        write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset+8*nglob) &
                            no_backward_acoustic_buffer(nglob*mod(no_backward_iframe+1,3)+1:nglob*(mod(no_backward_iframe+1,3)+1))

      endif
    endif ! end acoustic

    if (any_elastic) then

      !if (APPROXIMATE_HESS_KL) then
      !  offset = 4*2*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
      !else
      !  offset = 4*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
      !endif
     if(Full_Hessian_by_Wavefield_Stored) then
      if(SIMULATION_TYPE==1) then !---------------------------------------------------
         write(IOUT_UNDO_ATT) displ_elastic(:,:) 
         write(IOUT_UNDO_ATT+100) displ_elastic_m2(:,:)
         !-- for save accelaration for computing density kernels in Ha and Hb
         write(IOUT_UNDO_ATT+200) accel_elastic(:,:)
         write(IOUT_UNDO_ATT+300) accel_elastic_m2(:,:)
         close(IOUT_UNDO_ATT)
         close(IOUT_UNDO_ATT+100)
         close(IOUT_UNDO_ATT+200)
         close(IOUT_UNDO_ATT+300)
      else if(SIMULATION_TYPE==3) then !lucas, no accelaration needs to save here
         no_backward_displ_buffer(:,:) = displ_elastic_m2(:,:) 
         !write(6,*) 'it and no_backward_iframe',it,no_backward_iframe
         !write(6,*) 'no_backward_Nframes',no_backward_Nframes
         !write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset) no_backward_displ_buffer(:,:) ! lucas: for fortran2003
         write(IOUT_UNDO_ATT) no_backward_displ_buffer(:,:)
         close(IOUT_UNDO_ATT)
      endif!-----------------------------------------------------------------------------
     else !lucas, the old codes
       if(SIMULATION_TYPE==1) then !==========
        no_backward_displ_buffer(:,:) = displ_elastic(:,:) 
        write(IOUT_UNDO_ATT) no_backward_displ_buffer(:,:)
        write(IOUT_UNDO_ATT+200) accel_elastic(:,:)
        close(IOUT_UNDO_ATT)
        close(IOUT_UNDO_ATT+200)
       else if(SIMULATION_TYPE==3) then
         no_backward_displ_buffer(:,:) = displ_elastic(:,:) 
         write(IOUT_UNDO_ATT) no_backward_displ_buffer(:,:)
         close(IOUT_UNDO_ATT)
       endif !==================================
     endif

       
       if(write_field_at_frames) then!--------------------lucas, write few frames for plot wavefield (displ) by Matlab
        ! if(for_m1) then
        ! no_backward_displ_buffer(:,:) = displ_elastic(:,:) !for test
        ! endif
        ! if(.not. for_m1) then
        ! no_backward_displ_buffer(:,:) = displ_elastic_m2(:,:) !for test
        ! endif
       if(no_backward_iframe==000101 .or. no_backward_iframe==000301 .or. no_backward_iframe==000501  &
          .or. no_backward_iframe==000701 .or. no_backward_iframe==000901)then ! write one time frame for ploting,
          do ispec = 1, nspec
            do j = 1, 5
              do i = 1, 5
              iglob = ibool(i,j,ispec)
                 xx = coord(1,iglob)
                 zz = coord(2,iglob)
                 ! need to change when Full_Hessian_by_Wavefield_Stored=true
                 write(no_backward_iframe,'(4e15.4e4)') xx,zz,no_backward_displ_buffer(1,iglob),no_backward_displ_buffer(2,iglob)
              enddo
            enddo
          enddo
       endif
       endif ! end write_field_at_frames------------------

      !if (APPROXIMATE_HESS_KL) then
      !  no_backward_accel_buffer(:,:) = accel_elastic(:,:)
      !  write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset) no_backward_accel_buffer
      !endif
      !lucas---
      if(myrank==islice_selected_rec(1)) then
      iglob_source=ibool(2,2,ispec_selected_rec(1)) ! lucas, iglobe_source indiates the receiver
      if(mod(it,500)==0) then 
      write(6,*) 'print/check forward or adjoint fields at the first receiver location:'
      write(6,*) 'ispec_selected_source(1) = ',ispec_selected_rec(1)
      write(6,*) 'iglob_source = ',iglob_source
      write(6,*) 'displ_elastic(2,iglob_source)',displ_elastic(2,iglob_source)
      !write(6,*) 'displ_elastic_m2(2,iglob_source)',displ_elastic_m2(2,iglob_source)
      !write(6,*) 'b_displ_elastic(2,iglob_source)',b_displ_elastic(2,iglob_source)
      !write(6,*) 'b_displ_elastic_m2(2,iglob_source)',b_displ_elastic_m2(2,iglob_source)
      endif
     endif
     !lucas ---
    endif

 ! endif

  ! this operation will automatically synchronize the remaining I/O to do
  !if (it == NSTEP) close(IOUT_UNDO_ATT)
  !close(IOUT_UNDO_ATT)

  end subroutine save_forward_arrays_no_backward

