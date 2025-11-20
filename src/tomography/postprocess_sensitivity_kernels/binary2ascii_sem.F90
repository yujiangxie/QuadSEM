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
!Transfer models or kerenels from binary to ascii format
program binary2ascii_sem !lucas,

#ifdef USE_MPI
  use mpi
#endif

  use postprocess_par

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  integer, parameter :: NARGS = 3

  ! data must be of dimension: (NGLLX,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: dat
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat_store !lucas
  integer :: NGLOB_me, nspec_me
  ! MPI
  integer :: myrank,NPROC

  integer :: i,j,ier,ispec, iker
  character(len=MAX_STRING_LEN) :: arg(3)
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  character(len=MAX_STRING_LEN) :: prname

  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  !character(len=MAX_STRING_LEN) :: stations_filename,sources_filename !lucas
  integer :: nker

  character(len=MAX_STRING_LEN*2) :: ks_file
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: norm,max_old, min_old
  integer, dimension(:,:,:),allocatable :: ibool_me
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xstore_me, zstore_me
  real t1,t2
  ! MPI initialization
  ! lucas 1. initialize MPI
  call init_mpi()
  call world_size(NPROC)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running model_kernel_b2a_sem on",NPROC,"processors" !lucas
  call cpu_time(t1)

  ! parse command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
        print *, 'USAGE:  mpirun -np NPROC bin/xmodel_kernel_b2a_sem KERNEL_NAME INPUT_DIR OUPUT_DIR'
        call stop_the_code(' Please check command line arguments')
    endif
  endif

  ! synchronizes all processes
  call synchronize_all()
  !lucas 2. read command line input arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo

  kernel_names_comma_delimited = arg(1)
  input_dir= arg(2)
  output_dir = arg(3)
  
  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  allocate(norm(nker),max_old(nker),min_old(nker))

  ! user output
  if (myrank == 0) then
    print *,"command line arguments:"
    print *,"  input dir : ",trim(input_dir)
    print *,"  output dir: ",trim(output_dir)
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

  !lucas 4.
  ! frees arrays
  deallocate(dat,ibool_me)
 !---------------------------------------------------
 ! file output
 ! preconditioned kernel file name
   write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//'rhop_alpha_beta_kernel.dat'
   open(IOUT,file=trim(ks_file),status='unknown',iostat=ier)
   if (ier /= 0) call stop_the_code('Error opening precondtioned kernel file')
    do ispec = 1, nspec_me
      do j = 1, NGLLZ
        do i = 1, NGLLX
         write(IOUT,'(5e15.5e4)') xstore_me(i,j,ispec), zstore_me(i,j,ispec), &
                     dat_store(i,j,ispec,1),dat_store(i,j,ispec,2),dat_store(i,j,ispec,3)
        enddo
      enddo
    enddo
   close(IOUT)
   print *,'written: ',trim(ks_file) ! lucas, write for each processor
  !---------------------------------------------------


  ! frees memory
  deallocate(dat_store,xstore_me,zstore_me)

  ! synchronizes all processes
  call synchronize_all()
  !lucas 7. done!
#ifdef USE_MPI
  if (NPROC > 1) then

    ! the maximum value for the smoothed kernel
    norm(:) = max_old(:)
    call MPI_REDUCE(norm,max_old,nker,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    norm(:) = min_old(:)
    call MPI_REDUCE(norm,min_old,nker,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)


  endif
#endif

  do iker= 1, nker
    if (myrank == 0) then
      print *
      print *,' Min/Max for model or kernel = ',min_old(iker),max_old(iker), 'for ', trim(kernel_names(iker))


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

end program binary2ascii_sem

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
