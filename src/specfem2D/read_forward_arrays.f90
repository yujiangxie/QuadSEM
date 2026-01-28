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


  subroutine read_forward_arrays() !lucas, inlude CTD-SEM

! restores last time snapshot saved for backward/reconstruction of wavefields
! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
!          and adjoint sources will become more complicated
!          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields

  use constants, only: OUTPUT_FILES
  use specfem_par
  use specfem_par_gpu

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! acoustic medium
  if (any_acoustic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_acoustic**.bin')

    read(55) b_potential_acoustic
    read(55) b_potential_dot_acoustic
    read(55) b_potential_dot_dot_acoustic

    close(55)

    if (GPU_MODE) then
      ! transfers fields onto GPU
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic, &
                                          b_potential_dot_dot_acoustic, &
                                          Mesh_pointer)
    else
      ! free surface for an acoustic medium
      call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic)
    endif
  endif

  ! elastic medium
  if (any_elastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_elastic**.bin')

    read(55) b_displ_elastic
    read(55) b_veloc_elastic
    read(55) b_accel_elastic
    close(55)

    !SH (membrane) waves
    if (.not. P_SV) then
      ! only index array(1,:) contains SH wavefield
      b_displ_elastic(2,:) = 0._CUSTOM_REAL
      b_veloc_elastic(2,:) = 0._CUSTOM_REAL
      b_accel_elastic(2,:) = 0._CUSTOM_REAL
    endif
    if(CTD_SEM) then !lucas, read last state, CTD-SEM------start-----------------------
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic_m2_',myrank,'.bin'
    open(unit=550,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_elastic_m2_**.bin')

    read(550) b_displ_elastic_m2
    read(550) b_veloc_elastic_m2
    read(550) b_accel_elastic_m2
    close(550)

    !SH (membrane) waves
    if (.not. P_SV) then
      ! only index array(1,:) contains SH wavefield
      b_displ_elastic_m2(2,:) = 0._CUSTOM_REAL
      b_veloc_elastic_m2(2,:) = 0._CUSTOM_REAL
      b_accel_elastic_m2(2,:) = 0._CUSTOM_REAL
    endif
    endif!------------------CTD-SEM------end-------------------------   

    if (GPU_MODE) then
      ! prepares wavefields for transfering
      if (P_SV) then
        tmp_displ_2D(1,:) = b_displ_elastic(1,:)
        tmp_displ_2D(2,:) = b_displ_elastic(2,:)
        tmp_veloc_2D(1,:) = b_veloc_elastic(1,:)
        tmp_veloc_2D(2,:) = b_veloc_elastic(2,:)
        tmp_accel_2D(1,:) = b_accel_elastic(1,:)
        tmp_accel_2D(2,:) = b_accel_elastic(2,:)
      else
        ! SH waves
        tmp_displ_2D(1,:) = b_displ_elastic(1,:)
        tmp_displ_2D(2,:) = 0._CUSTOM_REAL
        tmp_veloc_2D(1,:) = b_veloc_elastic(1,:)
        tmp_veloc_2D(2,:) = 0._CUSTOM_REAL
        tmp_accel_2D(1,:) = b_accel_elastic(1,:)
        tmp_accel_2D(2,:) = 0._CUSTOM_REAL
      endif
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,Mesh_pointer)
    endif
  endif

  ! poroelastic medium
  if (any_poroelastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_poroelastic_s**.bin')

    read(55) b_displs_poroelastic
    read(55) b_velocs_poroelastic
    read(55) b_accels_poroelastic
    close(55)

    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_poroelastic_w**.bin')

    read(56) b_displw_poroelastic
    read(56) b_velocw_poroelastic
    read(56) b_accelw_poroelastic
    close(56)

    ! safety check
    if (GPU_MODE) then
      call stop_the_code('GPU_MODE error: sorry, reading lastframe from poroelastic simulation not implemented yet')
    endif
  endif

  end subroutine read_forward_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_undoatt() ! lucas added for m2 also

! reads in saved wavefields

  use constants, only: IIN_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,NGLLX,NGLLZ

  use specfem_par, only: myrank,iteration_on_subset,NSUBSET_ITERATIONS, &
    any_acoustic,any_elastic,ATTENUATION_VISCOACOUSTIC,ATTENUATION_VISCOELASTIC, &
    b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
    b_displ_elastic,b_veloc_elastic,b_accel_elastic,b_e1,b_e11,b_e13, &
    b_dux_dxl_old,b_duz_dzl_old,b_dux_dzl_plus_duz_dxl_old,b_e1_acous_sf,b_sum_forces_old, &
    b_displ_elastic_m2,b_veloc_elastic_m2,b_accel_elastic_m2,b_e1_m2,b_e11_m2,b_e13_m2, & !lucas, CTD-SEM
    b_dux_dxl_old_m2,b_duz_dzl_old_m2,b_dux_dzl_plus_duz_dxl_old_m2,CTD_SEM, & !lucas, CTD-SEM
    GPU_MODE,nspec_ATT_ac,nglob

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! current subset iteration
  iteration_on_subset_tmp = NSUBSET_ITERATIONS - iteration_on_subset + 1

!  write(6,*) 'lucas reading for time step: iteration_on_subset_tmp = ',iteration_on_subset_tmp
!  write(6,*) 'lucas reading for NSUBSET_ITERATIONS = ',NSUBSET_ITERATIONS

  ! reads in saved wavefield
  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
  ! opens corresponding snapshot file for reading
  open(unit=IIN_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for reading')

  if(CTD_SEM) then ! lucas add for reading m2 fields
  ! reads in saved wavefield
  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_m2_save_frame_at',iteration_on_subset_tmp,'.bin'
  ! opens corresponding snapshot file for reading
  open(unit=IIN_UNDO_ATT+1,file=trim(OUTPUT_FILES)//outputname, &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for reading')
  endif !---------------------------------------

  if (any_acoustic) then
    read(IIN_UNDO_ATT) b_potential_dot_dot_acoustic
    read(IIN_UNDO_ATT) b_potential_dot_acoustic
    read(IIN_UNDO_ATT) b_potential_acoustic
    if (GPU_MODE) call transfer_b_fields_ac_to_device(nglob,b_potential_acoustic,b_potential_dot_acoustic, &
                                                        b_potential_dot_dot_acoustic,Mesh_pointer)
    if (ATTENUATION_VISCOACOUSTIC) then
      read(IIN_UNDO_ATT) b_e1_acous_sf
      read(IIN_UNDO_ATT) b_sum_forces_old
      if (GPU_MODE) call transfer_viscoacoustic_b_var_to_device(NGLLX*NGLLZ*nspec_ATT_ac,b_e1_acous_sf, &
                                                                b_sum_forces_old,Mesh_pointer)
    endif
  endif

  if (any_elastic) then
    read(IIN_UNDO_ATT) b_accel_elastic
    read(IIN_UNDO_ATT) b_veloc_elastic
    read(IIN_UNDO_ATT) b_displ_elastic

    if (ATTENUATION_VISCOELASTIC) then
      read(IIN_UNDO_ATT) b_e1
      read(IIN_UNDO_ATT) b_e11
      read(IIN_UNDO_ATT) b_e13
      read(IIN_UNDO_ATT) b_dux_dxl_old
      read(IIN_UNDO_ATT) b_duz_dzl_old
      read(IIN_UNDO_ATT) b_dux_dzl_plus_duz_dxl_old
    endif
  endif
  close(IIN_UNDO_ATT)

  if(CTD_SEM) then ! lucas add for m2 ------
  if (any_elastic) then
    read(IIN_UNDO_ATT+1) b_accel_elastic_m2
    read(IIN_UNDO_ATT+1) b_veloc_elastic_m2
    read(IIN_UNDO_ATT+1) b_displ_elastic_m2

    if (ATTENUATION_VISCOELASTIC) then
      read(IIN_UNDO_ATT+1) b_e1_m2
      read(IIN_UNDO_ATT+1) b_e11_m2
      read(IIN_UNDO_ATT+1) b_e13_m2
      read(IIN_UNDO_ATT+1) b_dux_dxl_old_m2
      read(IIN_UNDO_ATT+1) b_duz_dzl_old_m2
      read(IIN_UNDO_ATT+1) b_dux_dzl_plus_duz_dxl_old_m2
    endif
  endif
  close(IIN_UNDO_ATT+1)
  endif !------

  end subroutine read_forward_arrays_undoatt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_no_backward()

  use constants, only: IIN_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,APPROXIMATE_HESS_KL,NDIM, CUSTOM_REAL !lucas: CUSTOM_REAL

  use specfem_par, only: myrank, any_acoustic,any_elastic,b_potential_acoustic, & 
    nglob,no_backward_acoustic_buffer,no_backward_displ_buffer, & 
    no_backward_iframe,no_backward_Nframes,GPU_MODE, &
    displ_elastic,b_displ_elastic,it,ibool,islice_selected_rec,ispec_selected_rec, b_displ_elastic_m2,&
    no_backward_displ_buffer_fwd_um2, no_backward_displ_buffer_adj_du_m, & !lucas
    no_backward_displ_buffer_fwd_du, no_backward_displ_buffer_adj_du_s, & !lucas
    Full_Hessian_by_Wavefield_Stored, no_backward_accel_buffer_fwd_du, & !lucas
    no_backward_accel_buffer_fwd_um1,no_backward_accel_buffer_fwd_um2,b_accel_elastic, & !lucas
    displ_elastic_m2, displ_elastic_m1 !lucas 
  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! local parameters
  integer :: ier1,ier2,buffer_num_async_IO,buffer_num_GPU_transfer ! lucas added ier2
  integer(KIND=8) :: offset
  character(len=MAX_STRING_LEN) :: outputname1, outputname2, outputname4, outputname5 !lucas, or use one filename
  !logical :: file_exists !lucas
  integer :: no_backward_reverse_iframe
  integer :: iglob_source_lucas  !lucas
          
  no_backward_displ_buffer_fwd_um2(:,:) = 0._CUSTOM_REAL  !lucas
  no_backward_displ_buffer_adj_du_m(:,:) = 0._CUSTOM_REAL  !lucas
  no_backward_displ_buffer_fwd_du(:,:) = 0._CUSTOM_REAL  !lucas
  no_backward_displ_buffer_adj_du_s(:,:) = 0._CUSTOM_REAL  !lucas

  ! lucas: below three lines copied from https://en.wikibooks.org/wiki/Fortran/Fortran_procedures_and_functions
  !intent(in) means that the variable value can enter, but not be changed
  !intent(out) means the variable is set inside the procedure and sent back to the main program with any initial values ignored.
  !intent(inout) means that the variable comes in with a value and leaves with a value (default). 

  ! EB EB June 2018 : in this routine, in order to overlap both disk = => RAM and RAM = => GPU transfers, we initiate the
  ! transfer of a wavefield two iterations before this wavefield is actually needed by the kernel computation.
  ! Two iterations before, the wavefield is read from the disk.
  ! One iteration before, this wavefield is transfered from the RAM to the GPU.
  ! In the text above, an iteration means NSTEP_BETWEEN_COMPUTE_KERNELS iterations of the timeloop.

  no_backward_iframe = no_backward_iframe + 1  ! lucas: this routine in manage_no_backward_reconstruction_io() of the time loop,  no_backward_iframe=0 set in prepare_timerun()
  ! lucas ----------------
  !if (it == 1) then
  !  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_No_backward_reconstruction_database.bin'
  !  ! opens corresponding file for reading
  !  INQUIRE(FILE=trim(OUTPUT_FILES)//outputname, EXIST=file_exists)
  !  write(6,*) 'lucas test file_exists =',file_exists 
  !  open(unit=IIN_UNDO_ATT,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname, &
  !    status='old',action='read',form='unformatted',access='stream',iostat=ier)
  ! write(6,*) 'ier and iostat ###############',ier
  ! if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_No_backward_reconstruction_database.bin for reading')
  !else
  !  wait(IIN_UNDO_ATT)
  !endif 
  ! lucas ----------------the no_backward_iframe=1 to NSTEP, no_backward_Nframes=NSTEP
    no_backward_reverse_iframe = no_backward_Nframes - no_backward_iframe + 2  ! lucas +2 for delta=10*dt, +1 for delta=1*dt

!   write(6,*) 'no_backward_Nframes, no_backward_iframe, no_backward_reverse_iframe=',no_backward_Nframes,&
!                                                             no_backward_iframe,no_backward_reverse_iframe
  
  if(Full_Hessian_by_Wavefield_Stored) then
      !lucas, read for u(m1) and u(m2) for computing Ha(du,u*), where du=u(m2)-u(m1) ----------used in m1
      write(outputname1,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m1_iframe_at',no_backward_reverse_iframe,'.bin'
      open(unit=IIN_UNDO_ATT,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname1, &
      status='old',action='read',form='unformatted',iostat=ier1)
      if (ier1 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_m1_iframe_at** for reading')
      !write(6,*) 'read proc***_disp_m1_iframe_at** successful'

      write(outputname2,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m2_iframe_at',no_backward_reverse_iframe,'.bin'
      open(unit=IIN_UNDO_ATT+1,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname2, &
      status='old',action='read',form='unformatted',iostat=ier2)
      if (ier2 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_m2_iframe_at** for reading')
      !write(6,*) 'read proc***_disp_m2_iframe_at** successful'

      !lucas read u(m1) and u*(m2) for computing Hb(u,du*), where du* = u*(m2) - u*(m1)----------used in m1
!      write(outputname3,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m2_iframe_at',no_backward_iframe,'.bin'  ! lucas: keep in mind that the u*(m) is used as its save
!      open(unit=IIN_UNDO_ATT+2,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname3, &
!      status='old',action='read',form='unformatted',iostat=ier2)
!      if (ier2 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m2_iframe_at** for reading')

      !for test
     ! write(outputname4,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m1_iframe_at',no_backward_iframe,'.bin'  ! lucas: keep in mind that the u*(m) is used as its save
     ! open(unit=IIN_UNDO_ATT+3,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname4, &
     ! status='old',action='read',form='unformatted',iostat=ier2)
     ! if (ier2 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m1_iframe_at** for reading')
     
      ! to store accelaration for computing the density of Ha, and density of Hb--
      write(outputname4,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_accel_m1_iframe_at',no_backward_reverse_iframe,'.bin'
      open(unit=IIN_UNDO_ATT+3,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname4, &
      status='old',action='read',form='unformatted',iostat=ier1)
      if (ier1 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_accel_m1_iframe_at** for reading')
      !write(6,*) 'read proc***_accel_m1_iframe_at** successful'

      write(outputname5,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_accel_m2_iframe_at',no_backward_reverse_iframe,'.bin'
      open(unit=IIN_UNDO_ATT+4,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname5, &
      status='old',action='read',form='unformatted',iostat=ier2)
      if (ier2 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_accel_m2_iframe_at** for reading')
      !write(6,*) 'read proc***_accel_m2_iframe_at** successful'
      

  else ! the old codes
      write(outputname1,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_disp_m1_iframe_at',no_backward_reverse_iframe,'.bin'
      open(unit=IIN_UNDO_ATT,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname1, &
      status='old',action='read',form='unformatted',iostat=ier1)
      if (ier1 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_disp_iframe_at** for reading')

    !  write(outputname1,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m2_iframe_at',no_backward_iframe,'.bin'  ! lucas: keep in mind that the u*(m) is used as its save
    !  open(unit=IIN_UNDO_ATT+1,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname1, &
    !  status='old',action='read',form='unformatted',iostat=ier2)
    !  if (ier2 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m2_iframe_at** for reading')

     
    !  write(outputname1,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_adj_disp_m1_iframe_at',no_backward_iframe,'.bin'  ! lucas: keep in mind that the u*(m) is used as its save
    !  open(unit=IIN_UNDO_ATT+2,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname1, &
    !  status='old',action='read',form='unformatted',iostat=ier2)
    !  if (ier2 /= 0 ) call exit_MPI(myrank,'Error opening file proc***_adj_disp_m2_iframe_at** for reading')
  endif
 
    !INQUIRE(FILE=trim(OUTPUT_FILES)//outputname1, EXIST=file_exists)
    !write(6,*) 'lucas test file_exists =',file_exists 
    !write(6,*) 'ier =',ier
    !write(6,*) ' no_backward_Nframes and no_backward_reverse_iframe = ', &
    !                           no_backward_Nframes, no_backward_reverse_iframe

  buffer_num_async_IO = mod(no_backward_iframe+2,3)
  buffer_num_GPU_transfer = mod(no_backward_iframe+1,3)

  if (any_acoustic) then

    ! offset is computed in two times to avoid integer overflow
    offset = 4*nglob
    offset = offset*(no_backward_Nframes - no_backward_iframe ) + 1
    if (no_backward_iframe <= no_backward_Nframes) &
      read(IIN_UNDO_ATT,asynchronous='yes',pos=offset) &
           no_backward_acoustic_buffer(nglob*buffer_num_async_IO+1:nglob*(buffer_num_async_IO+1))

    if (no_backward_iframe /= 1) then
      if (GPU_MODE) then
        ! this function ensures the previous async transfer is finished, and
        ! launches the transfer of the next wavefield
        call transfer_async_pot_ac_to_device(no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1),Mesh_pointer)
      else
        ! we get the wavefield from the previous iteration, because this RAM = =>
        ! RAM copy is blocking
        b_potential_acoustic(:) = no_backward_acoustic_buffer(nglob*mod(no_backward_iframe,3)+1: &
                                                                nglob*(mod(no_backward_iframe,3)+1))
      endif
    endif

  endif ! any_acoustic

  if (any_elastic) then
     
      if(Full_Hessian_by_Wavefield_Stored) then
        read(IIN_UNDO_ATT) no_backward_displ_buffer ! lucas: u(m1)
        read(IIN_UNDO_ATT+1) no_backward_displ_buffer_fwd_um2 ! lucas read u(m2) for computing du
        !for du
        no_backward_displ_buffer_fwd_du = no_backward_displ_buffer_fwd_um2 - no_backward_displ_buffer ! du=u(m2)-u(m1)
        !for du*
        no_backward_displ_buffer_adj_du_s = displ_elastic_m1 - displ_elastic   ! du*=u*_s(m1) - u*(m1)
        no_backward_displ_buffer_adj_du_m = displ_elastic_m2 - displ_elastic   ! du*=u*_m(m2) - u*(m1)

       ! read accel for density kernels Ha and Hb
        read(IIN_UNDO_ATT+3) no_backward_accel_buffer_fwd_um1 ! lucas read accel(m1) 
        read(IIN_UNDO_ATT+4) no_backward_accel_buffer_fwd_um2 ! lucas read accel*(m2) 
        !lucas, for Ha in density kernels
        no_backward_accel_buffer_fwd_du=no_backward_accel_buffer_fwd_um2-no_backward_accel_buffer_fwd_um1 
        !lucas, for Hb in density kernels
        b_accel_elastic=no_backward_accel_buffer_fwd_um1 
        !-----
        b_displ_elastic_m2=no_backward_displ_buffer_fwd_um2           
        b_displ_elastic=   no_backward_displ_buffer 
   
      else !lucas, old codes
        read(IIN_UNDO_ATT) no_backward_displ_buffer ! lucas: u(m1)
         b_displ_elastic= no_backward_displ_buffer !u(m1) = read in
      endif

      !=============================user ouput for check==============================
      if(Full_Hessian_by_Wavefield_Stored) then !no SIMULATION_TYPE==1 in this code
      ! 3.1 lucas: test for du* 
      if(myrank==islice_selected_rec(1)) then
        iglob_source_lucas=ibool(2,2,ispec_selected_rec(1)) ! lucas: find the iglobe of element(2,2) at the receiver.
        if(mod(it,500)==0) then
          write(6,*) 'only valid in du* check at the first receiver location, there are du*, u*(m), u*(m+dm):' 
          write(6,*) 'ispec_selected_source(1) = ',ispec_selected_rec(1)
          write(6,*) 'iglob_source_lucas = ',iglob_source_lucas
          write(6,*) 'no_backward_displ_buffer_adj_du_m(2,iglob_source_lucas) = du*_m = ', &
                      no_backward_displ_buffer_adj_du_m(2,iglob_source_lucas) ! lucas: du*_m
          write(6,*) 'displ_elastic(2,iglob_source_lucas) = u*(m) = ',displ_elastic(2,iglob_source_lucas) !lucas: u*(m)
          write(6,*) 'displ_elastic_m2(2,iglob_source_lucas)= u*(m+dm) = ', & 
                      displ_elastic_m2(2,iglob_source_lucas) ! lucas: u*(m+dm)
          write(6,*) '=============================================='
        endif
      endif
      endif
     !------------------------------
     if(Full_Hessian_by_Wavefield_Stored) then ! for model m
     ! 3.2 lucas: test for du , the tested point should be the receiver or the source
     if(myrank==islice_selected_rec(1)) then
       iglob_source_lucas=ibool(2,2,ispec_selected_rec(1))
       if(mod(it,500)==0) then 
       write(6,*) 'only valid in du check at the first receiver location, there are du, u(m), u(m+dm):'
       write(6,*) 'ispec_selected_source(1) = ',ispec_selected_rec(1)
       write(6,*) 'iglob_source_lucas = ',iglob_source_lucas
       write(6,*) 'no_backward_displ_buffer_fwd_du(2,iglob_source_lucas) = du = ', & 
                   no_backward_displ_buffer_fwd_du(2,iglob_source_lucas) !lucas: du
       write(6,*) 'no_backward_displ_buffer(2,iglob_source_lucas) = u(m) =', & 
                      no_backward_displ_buffer(2,iglob_source_lucas) ! lucas: u(m)
       write(6,*) 'no_backward_displ_buffer_fwd_um2(2,iglob_source_lucas) = u(m+dm) =', & 
                      no_backward_displ_buffer_fwd_um2(2,iglob_source_lucas) ! lucas: u(m+dm)
       endif
     endif
     endif
  endif ! end any_elastic
  
  if(Full_Hessian_by_Wavefield_Stored) then
    close(IIN_UNDO_ATT)
    close(IIN_UNDO_ATT+1)
    close(IIN_UNDO_ATT+2)
    close(IIN_UNDO_ATT+3)
    close(IIN_UNDO_ATT+4)
  else
    close(IIN_UNDO_ATT)
  endif
  ! lucas ----------------------------------------------------------------------------------------------------

   ! start old version here---------------------------------------------------------
   ! if (any_elastic) then
   !    b_displ_elastic(:,:) = no_backward_displ_buffer(:,:)
   !    if (APPROXIMATE_HESS_KL) b_accel_elastic(:,:) = no_backward_accel_buffer(:,:)
   !    if (APPROXIMATE_HESS_KL) then
   !     offset = 4*2*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
   !    else
   !     offset = 4*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
   !    endif
   !    read(IIN_UNDO_ATT,asynchronous='yes',pos=offset) no_backward_displ_buffer(:,:)
   !    if (APPROXIMATE_HESS_KL) read(IIN_UNDO_ATT,asynchronous='yes',pos=offset+8+4*nglob) no_backward_accel_buffer(:,:)
   ! endif
   ! if (it == NSTEP) close(IIN_UNDO_ATT)
   ! end old version here------------------------------------------------------------

  end subroutine read_forward_arrays_no_backward

