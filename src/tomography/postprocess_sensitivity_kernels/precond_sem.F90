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
!lucas: Yujiang Xie, changed from program smooth_sem
!lucas:USAGE:  mpirun -np NPROC bin/xsmooth_sem a_s b_s a_r b_r KERNEL_NAME INPUT_DIR OUPUT_DIR binary_kernel(F or T)
!lucas: KERNEL_NAME only support for rhop_kernel, alpha_kernel, beta_kernel, but can be easy developed for kappa_kernel, etc.
!lucas: a_s, b_s for sources, and a_r, b_r for receivers. e.g., precond_s =1/(1+a_s*exp(-b_s*r)); precond_r=1/(1+a_r*exp(-b_r*r^2))
program precond_sem !lucas, I changed this to precond for sources and receivers

#ifdef USE_MPI
  use mpi
#endif

  use postprocess_par

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  integer, parameter :: NARGS = 8

  ! data must be of dimension: (NGLLX,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: dat
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: dist_sta, dist_src, precond_receiver, precond_source !lucas
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat_store,dat_precond !lucas
  integer :: NGLOB_me, nspec_me
  ! MPI
  integer :: myrank,NPROC

  integer :: i,j,ier,ispec, iker, iglob
  integer :: irec, nrec, isrc, nsrc !lucas, for receivers and sources

  character(len=MAX_STRING_LEN) :: arg(8)
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  character(len=MAX_STRING_LEN) :: prname

  logical :: binary_kernel

  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: stations_filename,sources_filename !lucas
  integer :: nker

  ! smoothing parameters
  character(len=MAX_STRING_LEN*2) :: ks_file
  !real(kind=CUSTOM_REAL) :: sigma_h, sigma_v
  !real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2_inv, sigma_h3_sq, sigma_v,sigma_v2_inv, sigma_v3_sq
  !real(kind=CUSTOM_REAL) :: norm_h, norm_v
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: norm,max_old,max_new, min_old, min_new
  real(kind=CUSTOM_REAL) :: a_r, b_r, a_s, b_s !lucas, for preconditioner operators
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: exp_val, factor, wgll_sq
  !double precision, dimension(NGLLX) :: wxgll
  !double precision, dimension(NGLLX) :: xigll
  integer, dimension(:,:,:),allocatable :: ibool_me
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xstore_me, zstore_me
  !real(kind=CUSTOM_REAL) :: dist_h,dist_v
  !real(kind=CUSTOM_REAL) :: element_size
  !real(kind=CUSTOM_REAL), PARAMETER :: PI = 3.1415927
  real t1,t2
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x_sta,z_sta
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x_src,z_src
  ! MPI initialization
  ! lucas 1. initialize MPI
  call init_mpi()
  call world_size(NPROC)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running XPRECOND_SEM on",NPROC,"processors" !lucas
  call cpu_time(t1)

  ! parse command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
        print *, 'USAGE:  mpirun -np NPROC bin/xprecond_sem a_s b_s a_r b_r KERNEL_NAME INPUT_DIR OUPUT_DIR binary_kernel(F or T),'
        call stop_the_code(' Please check command line arguments')
    endif
  endif

  ! synchronizes all processes
  call synchronize_all()
  !lucas 2. read command line input arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo

  read(arg(1),*) a_s
  read(arg(2),*) b_s
  read(arg(3),*) a_r
  read(arg(4),*) b_r
  kernel_names_comma_delimited = arg(5)
  input_dir= arg(6)
  output_dir = arg(7)
  read(arg(8),*) binary_kernel

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  allocate(norm(nker),max_new(nker),max_old(nker),min_new(nker),min_old(nker))
 !if (GPU_MODE) call initialize_cuda_device(myrank,ncuda_devices)


  ! user output
  if (myrank == 0) then
    print *,"command line arguments:"
    print *,"  preconditioner parameters: a_s , b_s for sources: ",a_s,b_s
    print *,"  preconditioner parameters: a_r , b_r for receiver: ",a_r,b_r
    ! scalelength: approximately S ~ sigma * sqrt(8.0) for a Gaussian smoothing
    !print *,"  smoothing scalelengths horizontal, vertical: ",sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *,"  input dir : ",trim(input_dir)
    print *,"  output dir: ",trim(output_dir)
    print *,"  output binary kernel or not, binary_kernel = ", binary_kernel
    print *
  endif

  !lucas 3. read input files for considered Gll points, ibool, _x, _z, kernels
  write(prname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'/proc',myrank,'_NSPEC_ibool.bin'
  open(IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error opening _NSPEC_IBOOL file')
  endif
  read(IIN) nspec_me
  allocate(ibool_me(NGLLX,NGLLZ,nspec_me))
  read(IIN) ibool_me
  close(IIN)
  nglob_me = maxval(ibool_me(:,:,:))
  allocate(xstore_me(NGLLX,NGLLZ,NSPEC_me),zstore_me(NGLLX,NGLLZ,NSPEC_me),stat=ier)

  write(prname, '(a,i6.6,a)') trim(IN_DATA_FILES)//'/proc',myrank,'_x.bin'
    ! gets the coordinate x of the points located in my slice
  open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error reading neighbors external mesh file')
  endif
    ! global point arrays
  read(IIN) xstore_me
  close(IIN)

  write(prname, '(a,i6.6,a)') trim(IN_DATA_FILES)//'/proc',myrank,'_z.bin'
    ! gets the coordinate z of the points located in my slice
  open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error reading neighbors external mesh file')
  endif
    ! global point arrays
  read(IIN) zstore_me
  close(IIN)

  ! read sources and stations locations ===========================
  stations_filename=trim(IN_DATA_FILES)//'STATIONS_all'
  sources_filename=trim(IN_DATA_FILES)//'SOURCES_all'
  nrec=39
  nsrc=16
  allocate(x_sta(nrec),z_sta(nrec))
  allocate(x_src(nsrc),z_src(nsrc))

  open(unit= 1,file=trim(stations_filename),status='old',action='read')
  do irec = 1,nrec
  read(1,*) x_sta(irec),z_sta(irec)
  enddo
  close(1)

  open(unit= 2,file=trim(sources_filename),status='old',action='read')
  do isrc = 1,nsrc
  read(2,*) x_src(isrc),z_src(isrc)
  enddo
  close(2)

  allocate(dat(NGLLX,NGLLZ,NSPEC_me),dat_store(NGLLX,NGLLZ,NSPEC_me,nker),stat=ier) !lucas: nker could be vp_kernel, vs_kernel, etc
  if (ier /= 0) call stop_the_code('Error allocating dat array')
  do iker= 1, nker
      ! data file
    write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(kernel_names(iker))//'.bin'
    open(unit = IIN,file = trim(prname),status='old',action='read',form ='unformatted',iostat=ier)
    if (ier /= 0) then
        print *,'Error opening data file: ',trim(prname)
        call stop_the_code('Error opening data file')
    endif
    read(IIN) dat ! (NGLLX, NGLLZ, nspec)
    close(IIN)

    dat_store(:,:,:,iker) = dat(:,:,:)

    
    max_old(iker) = maxval(abs(dat(:,:,:)))
    min_old(iker) = minval(abs(dat(:,:,:)))
    
  enddo

  ! synchronizes all processes
  call synchronize_all()

  ! GPU setup
  !if (GPU_MODE) then
  !  call prepare_arrays_GPU(Container,xstore_me,zstore_me, &
  !                          sigma_h2_inv,sigma_v2_inv,sigma_h3_sq,sigma_v3_sq,nspec_me,nker,wgll_sq)
  !  ! synchronizes all processes
  !  call synchronize_all()
  !endif

  !lucas 4. initialize dat_precond and the parameters for the preconditioner operators
  allocate(dat_precond(NGLLX,NGLLZ,NSPEC_me,nker),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array dat_precond')
  dat_precond(:,:,:,:) = 0.0_CUSTOM_REAL
  
  ! lucas 5. to get dat_precond for receivers and sources
  !if (GPU_MODE) then
  !  call compute_smooth(Container,jacobian,xstore_other,zstore_other,dat_store,nspec_other)
  !endif
  
  allocate(dist_sta(NGLLX,NGLLZ,NSPEC_me),precond_receiver(NGLLX,NGLLZ,NSPEC_me),stat=ier)
  !lucas, preconditioner_station applies to all Gll points of current processor
  do irec = 1,nrec !------------------loop for all receivers
    do ispec = 1, nspec_me
       !----
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool_me(i,j,ispec)
            dist_sta(i,j,ispec) = sqrt((x_sta(irec)-xstore_me(i,j,ispec))**2 + (z_sta(irec)-zstore_me(i,j,ispec))**2) !lucas distant to stations, i.e., r in the formula
            precond_receiver(i,j,ispec)=1/(1 + a_r*exp(-b_r*dist_sta(i,j,ispec)*dist_sta(i,j,ispec)))
               do iker= 1, nker
                 dat_precond(i,j,ispec,iker)=precond_receiver(i,j,ispec)*dat_store(i,j,ispec,iker)
               enddo

          enddo
        enddo
       !----
    enddo ! ispec
    
    ! loop over each station using previous station's preconditioned kernel
   do ispec = 1, nspec_me
    do j = 1,NGLLZ
      do i = 1,NGLLX
        do iker= 1, nker
           dat_store(i,j,ispec,iker)=dat_precond(i,j,ispec,iker)
        enddo
      enddo
    enddo
   enddo

  enddo ! nrec -------------------
  !========================================================================================================
  allocate(dist_src(NGLLX,NGLLZ,NSPEC_me),precond_source(NGLLX,NGLLZ,NSPEC_me),stat=ier)
  !lucas, preconditioner_stations applies to all Gll points of current processor
  do isrc = 1,nsrc !------------------loop for all stations
    do ispec = 1, nspec_me
       !----
       do j = 1,NGLLZ
         do i = 1,NGLLX
          iglob = ibool_me(i,j,ispec)
          dist_src(i,j,ispec) = sqrt((x_src(isrc)-xstore_me(i,j,ispec))**2 + (z_src(isrc)-zstore_me(i,j,ispec))**2) !lucas distant to sources, i.e., r in the formula
          precond_source(i,j,ispec)=1/(1 + a_s*exp(-b_s*dist_src(i,j,ispec)))
            do iker= 1, nker
              dat_precond(i,j,ispec,iker)=precond_source(i,j,ispec)*dat_store(i,j,ispec,iker)
            enddo

         enddo
       enddo
       !----
    enddo ! ispec

    ! loop over each station using previous station's preconditioned result
   do ispec = 1, nspec_me
    do j = 1,NGLLZ
      do i = 1,NGLLX
        do iker= 1, nker
           dat_store(i,j,ispec,iker)=dat_precond(i,j,ispec,iker)
        enddo
      enddo
    enddo
   enddo

  enddo ! nsrc -------------------

  !lucas 6.
  ! frees arrays
  deallocate(dat,dat_store,dist_sta,dist_src,precond_receiver,precond_source,ibool_me)

  do iker =1, nker
   max_new(iker) = maxval(abs(dat_precond(:,:,:,iker)))
   min_new(iker) = minval(abs(dat_precond(:,:,:,iker)))
  enddo

if(binary_kernel) then!lucas ++++++++++
  do iker= 1, nker
    ! file output
    ! smoothed kernel file name
    write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_names(iker))//'_precond.bin'
    open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening smoothed kernel file')
    write(IOUT) dat_precond(:,:,:,iker)
    close(IOUT)
    print *,'written: ',trim(ks_file) ! lucas, write for each processor
  enddo
else
 ! lucas 6.1 save for ascii
 !---------------------------------------------------
 ! file output
 ! preconditioned kernel file name
   write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//'rhop_alpha_beta_kernel'//'_precond.dat'
   open(IOUT,file=trim(ks_file),status='unknown',iostat=ier)
   if (ier /= 0) call stop_the_code('Error opening precondtioned kernel file')
    do ispec = 1, nspec_me
      do j = 1, NGLLZ
        do i = 1, NGLLX
         write(IOUT,'(5e15.5e4)') xstore_me(i,j,ispec), zstore_me(i,j,ispec), &
                     dat_precond(i,j,ispec,1),dat_precond(i,j,ispec,2),dat_precond(i,j,ispec,3)
        enddo
      enddo
    enddo
   close(IOUT)
   print *,'written: ',trim(ks_file) ! lucas, write for each processor
  !---------------------------------------------------
endif!lucas ++++++++++


  ! frees memory
  deallocate(dat_precond)
  deallocate(xstore_me,zstore_me)

  ! synchronizes all processes
  call synchronize_all()
  !lucas 7. done!
#ifdef USE_MPI
  if (NPROC > 1) then

    ! the maximum value for the smoothed kernel
    norm(:) = max_old(:)
    call MPI_REDUCE(norm,max_old,nker,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    norm(:) = max_new(:)
    call MPI_REDUCE(norm,max_new,nker,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    norm(:) = min_old(:)
    call MPI_REDUCE(norm,min_old,nker,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

    norm(:) = min_new(:)
    call MPI_REDUCE(norm,min_new,nker,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  endif
#endif

  do iker= 1, nker
    if (myrank == 0) then
      print *
      print *,' Min / Max data value before preconditioner = ', min_old(iker), max_old(iker), 'for ', trim(kernel_names(iker))
      print *,' Min / Max data value after preconditioner  = ', min_new(iker), max_new(iker), 'for ', trim(kernel_names(iker))

    endif
  enddo

  call cpu_time(t2)

  !if (GPU_Mode) then
  !  print *,'Computation time with GPU:',t2-t1
  !else
    print *,'Computation time with CPU:',t2-t1
  !endif

  if (myrank == 0) close(IIN)

  ! MPI finish
  call finalize_mpi()

end program precond_sem

!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,z0,sigma_h2_inv,sigma_v2_inv,exp_val,xx_elem,zz_elem)

  use constants
  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ),intent(in) :: xx_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,z0,sigma_h2_inv,sigma_v2_inv


  ! local parameters
  integer :: ii,jj
  real(kind=CUSTOM_REAL) :: dist_h,dist_v


  do jj = 1, NGLLZ
    do ii = 1, NGLLX
      ! gets vertical and horizontal distance
      call get_distance_square_vec(dist_h,dist_v,x0,z0,xx_elem(ii,jj),zz_elem(ii,jj))
       ! Gaussian function
      exp_val(ii,jj) = exp(- sigma_h2_inv*dist_h - sigma_v2_inv*dist_v)
    enddo
  enddo

  end subroutine smoothing_weights_vec


!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_square_vec(dist_h,dist_v,x0,z0,x1,z1)

! returns square lengths as distances in radial and horizontal direction
! only for flat earth with z in vertical direction

  use constants
  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,z0,x1,z1

  ! vertical distance
  dist_v =  (z0-z1)**2

  ! horizontal distance
  dist_h =  (x0-x1)**2

  end subroutine get_distance_square_vec
