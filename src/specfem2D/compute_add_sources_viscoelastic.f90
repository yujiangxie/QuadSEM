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

! for viscoelastic solver

  subroutine compute_add_sources_viscoelastic(accel_elastic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: myrank,P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + sourcearrays(i_source,2,i,j) * stf_used
              
              !lucas debug
              if(it==100 .and. i==5 .and. j==5) then
              print *, '*********########********* at compute_add_sources_viscoelastic()'
              print *, '*********########*********sourcearrays(i_source,1,5,5)=',sourcearrays(i_source,1,5,5) 
              print *, '*********########*********sourcearrays(i_source,2,5,5)=',sourcearrays(i_source,2,5,5)
              print *, '*********########*********accel_elastic(1,iglob)=',accel_elastic(1,iglob)
              print *, '*********########*********accel_elastic(2,iglob)=',accel_elastic(2,iglob)
              print *, '*********########*********stf_used=',stf_used
              endif

            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used
            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic

! lucas add it for backward simulation test, not used
! subroutine compute_add_sources_viscoelastic_backward(b_accel_elastic,it,i_stage) 

!  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

!  use specfem_par, only: myrank,P_SV,ispec_is_elastic,nglob_elastic, &
!                         NSOURCES,source_time_function, &
!                         islice_selected_source,ispec_selected_source,sourcearrays, &
!                         ibool
!  implicit none

!  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: b_accel_elastic
!  integer :: it, i_stage

!  !local variable
!  integer :: i_source,i,j,iglob,ispec
!  real(kind=CUSTOM_REAL) :: stf_used

  ! --- add the source
!  do i_source = 1,NSOURCES

    ! if this processor core carries the source
!    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
!      ispec = ispec_selected_source(i_source)

      ! source element is elastic
!      if (ispec_is_elastic(ispec)) then

        ! source time function
!        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
!        if (P_SV) then
          ! P-SV calculation
!          do j = 1,NGLLZ
!            do i = 1,NGLLX
!              iglob = ibool(i,j,ispec)
!              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used
!              b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) + sourcearrays(i_source,2,i,j) * stf_used
              
!              !lucas debug
!              if(it==100 .and. i==5 .and. j==5) then
!              print *, '*********########********* at compute_add_sources_viscoelastic()'
!              print *, '*********########*********sourcearrays(i_source,1,5,5)=',sourcearrays(i_source,1,5,5) 
!              print *, '*********########*********sourcearrays(i_source,2,5,5)=',sourcearrays(i_source,2,5,5)
!              print *, '*********########*********b_accel_elastic(1,iglob)=',b_accel_elastic(1,iglob)
!              print *, '*********########*********b_accel_elastic(2,iglob)=',b_accel_elastic(2,iglob)
!              print *, '*********########*********stf_used=',stf_used
!              endif

!            enddo!
!          enddo
!        else
          ! SH (membrane) calculation
!          do j = 1,NGLLZ
!            do i = 1,NGLLX
!              iglob = ibool(i,j,ispec)
!              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used
!            enddo
!          enddo
!        endif

!      endif ! source element is elastic
!    endif ! if this processor core carries the source
!  enddo ! do i_source= 1,NSOURCES

!  end subroutine compute_add_sources_viscoelastic_backward





! lucas, CTD-SEM--------m2-----------------------------------

subroutine compute_add_sources_viscoelastic_m2(accel_elastic_m2,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: myrank,P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic_m2
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              !lucas, CTD-SEM ---------------------
               accel_elastic_m2(1,iglob) = accel_elastic_m2(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used
               accel_elastic_m2(2,iglob) = accel_elastic_m2(2,iglob) + sourcearrays(i_source,2,i,j) * stf_used
             
              !lucas debug
              if(it==100 .and. i==5 .and. j==5) then
              print *, '*********########********* at compute_add_sources_viscoelastic() in CTD-SEM'
              print *, '*********########*********sourcearrays(i_source,1,5,5)=',sourcearrays(i_source,1,5,5) 
              print *, '*********########*********sourcearrays(i_source,2,5,5)=',sourcearrays(i_source,2,5,5)
              print *, '*********########*********accel_elastic_m2(1,iglob)=',accel_elastic_m2(1,iglob)
              print *, '*********########*********accel_elastic_m2(2,iglob)=',accel_elastic_m2(2,iglob)
              print *, '*********########*********stf_used=',stf_used
              endif

            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
             
              accel_elastic_m2(1,iglob) = accel_elastic_m2(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used  ! lucas, CTD-SEM
             
            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic_m2


! lucas, CTD-SEM--------m1-----------------------------------

subroutine compute_add_sources_viscoelastic_m1(accel_elastic_m1,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: myrank,P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic_m1
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              !lucas, CTD-SEM ---------------------
               accel_elastic_m1(1,iglob) = accel_elastic_m1(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used
               accel_elastic_m1(2,iglob) = accel_elastic_m1(2,iglob) + sourcearrays(i_source,2,i,j) * stf_used
             
              !lucas debug
              if(it==100 .and. i==5 .and. j==5) then
              print *, '*********########********* at compute_add_sources_viscoelastic() in CTD-SEM'
              print *, '*********########*********sourcearrays(i_source,1,5,5)=',sourcearrays(i_source,1,5,5) 
              print *, '*********########*********sourcearrays(i_source,2,5,5)=',sourcearrays(i_source,2,5,5)
              print *, '*********########*********accel_elastic_m1(1,iglob)=',accel_elastic_m1(1,iglob)
              print *, '*********########*********accel_elastic_m1(2,iglob)=',accel_elastic_m1(2,iglob)
              print *, '*********########*********stf_used=',stf_used
              endif

            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
             
              accel_elastic_m1(1,iglob) = accel_elastic_m1(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used  ! lucas, CTD-SEM
             
            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic_m1


!
!=====================================================================
!

  subroutine compute_add_sources_viscoelastic_moving_source(accel_elastic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,SOURCE_IS_MOVING,TINYVAL,NGLJ,IMAIN

  use specfem_par, only: P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool,coord,nspec,nglob,xigll,zigll,z_source,NPROC, & !These 3 lines are for moving src
                         xi_source,gamma_source,coorg,knods,ngnod,npgeo,iglob_source,x_source,z_source,deltat,t0,myrank, &
                         time_stepping_scheme,hxis_store,hgammas_store,tshift_src,source_type,ispec_is_acoustic, &
                         hxis,hpxis,hgammas,hpgammas,anglesource,ispec_is_poroelastic,Mxx,Mxz,Mzz,gammax,gammaz,xix,xiz, &
                         AXISYM,xiglj,is_on_the_axis,initialfield
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used
  double precision :: hlagrange
  double precision :: xminSource,vSource,timeval,t_used
  ! single source array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  ! checks if anything to do
  if (.not. SOURCE_IS_MOVING) return

  !xminSource = -5000.0d0 !m
  !vSource = 2150.0d0 !1425.0d0 !1250.0 !m/s
  xminSource = -60.0d0 !m
  vSource = 60.0d0 !m/s

  if (time_stepping_scheme == 1) then
    ! Newmark
    timeval = (it-1)*deltat
  else
    call exit_MPI(myrank,'Not implemented!')
  endif

  ! moves and re-locates sources along x-axis

  do i_source = 1,NSOURCES
    if (abs(source_time_function(i_source,it,i_stage)) > TINYVAL) then
      t_used = (timeval-t0-tshift_src(i_source))

      x_source(i_source) = xminSource + vSource*t_used !timeval?

      ! collocated force source
      if (source_type(i_source) == 1) then
        call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           x_source(i_source),z_source(i_source), &
                           ispec_selected_source(i_source),islice_selected_source(i_source), &
                           NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                           iglob_source(i_source),.true.)

      else if (source_type(i_source) == 2) then
        ! moment-tensor source
        call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           x_source(i_source),z_source(i_source), &
                           ispec_selected_source(i_source),islice_selected_source(i_source), &
                           NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                           iglob_source(i_source),.false.)

      else if (.not. initialfield) then

        call exit_MPI(myrank,'incorrect source type')

      endif

      ispec = ispec_selected_source(i_source)
      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! Lagrange interpolators
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            call lagrange_any(xi_source(i_source),NGLJ,xiglj,hxis,hpxis)
            !do j = 1,NGLJ ! ABAB same result with that loop, this is good
            !  hxis(j) = hglj(j-1,xi_source(i_source),xiglj,NGLJ)
            !enddo
          else
            call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
          endif
        else
          call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
        endif
        call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)

        if (mod(it,10000) == 0) then
            !  write(IMAIN,*) "myrank:",myrank
            ! user output
            if (myrank == islice_selected_source(i_source)) then
              iglob = ibool(2,2,ispec_selected_source(i_source))
              !write(IMAIN,*) 'xcoord: ',coord(1,iglob)
              write(IMAIN,*) 'Problem... it??: ',it,'xcoord: ',coord(1,iglob)," iglob",iglob
              !'source carried by proc',myrank,"  source x:",x_source(i_source)," ispec:",ispec_selected_source(i_source)

              !call flush_IMAIN()
            endif

        endif

        ! stores Lagrangians for source
        hxis_store(i_source,:) = hxis(:)
        hgammas_store(i_source,:) = hgammas(:)

        sourcearray(:,:,:) = 0._CUSTOM_REAL

        ! computes source arrays
        select case (source_type(i_source))
        case (1)
          ! collocated force source
          do j = 1,NGLLZ
            do i = 1,NGLLX
              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)

              ! source element is acoustic
              if (ispec_is_acoustic(ispec)) then
                sourcearray(:,i,j) = real(hlagrange,kind=CUSTOM_REAL)
              endif

              ! source element is elastic
              if (ispec_is_elastic(ispec)) then
                if (P_SV) then
                  ! P_SV case
                  sourcearray(1,i,j) = real(- sin(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                  sourcearray(2,i,j) = real(cos(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                else
                  ! SH case (membrane)
                  sourcearray(:,i,j) = real(hlagrange,kind=CUSTOM_REAL)
                endif
              endif

              ! source element is poroelastic
              if (ispec_is_poroelastic(ispec)) then
                sourcearray(1,i,j) = real(- sin(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                sourcearray(2,i,j) = real(cos(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
              endif

            enddo
          enddo

        case (2)
          ! moment-tensor source
          call compute_arrays_source(ispec,xi_source(i_source),gamma_source(i_source),sourcearray, &
                                     Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)
          ! checks source
          if (ispec_is_acoustic(ispec)) then
            call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
          endif

          ! checks wave type
          if (ispec_is_elastic(ispec)) then
            if (.not. P_SV ) call exit_MPI(myrank,'cannot have moment tensor source in SH (membrane) waves calculation')
          endif

        end select

        ! stores sourcearray for all sources
        sourcearrays(i_source,:,:,:) = sourcearray(:,:,:)

      endif
    endif
  enddo

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                                       sourcearrays(i_source,1,i,j) * stf_used
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + &
                                       sourcearrays(i_source,2,i,j) * stf_used
            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                                       sourcearrays(i_source,1,i,j) * stf_used

              ! daniel debug source contribution
              !if (iglob == 37905) &
              !write(1234,*) it, dble(sourcearrays(i_source,1,i,j) * source_time_function(i_source,it,i_stage)), &
              !              accel_elastic(1,iglob),source_time_function(i_source,it,i_stage),sourcearrays(i_source,1,i,j)


            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic_moving_source

!
!=====================================================================
!

! for viscoelastic solver for adjoint propagation wave field

  subroutine compute_add_sources_viscoelastic_adjoint() 

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM !lucas added NDIM
  

  use specfem_par, only: P_SV,accel_elastic, ispec_is_elastic,NSTEP,it, & 
                         nrecloc,ispec_selected_rec_loc,ibool, & !myrank, & ! lucas added myrank 
                         source_adjoint,xir_store_loc,gammar_store_loc, anglesource
                         !lucas, xir_store_loc is the lagrange interpolant of receiver, see setup_source_receivers.f90, !lucas add anglesource and sourcearrays
                        
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec,i_source ! lucas added i_source=1
  integer :: it_tmp !xx,yy
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearrays_lucas

  ! time step index
  it_tmp = NSTEP - it + 1
              i_source=1
             
             ! if(it_tmp==7200-100) then
             ! do yy=1,NGLLZ
             ! do xx=1,NGLLX
             ! print *, '*********########********* before the do'
             ! print *, '*********########********* myrank=',myrank
             ! print *, '*********########*********sourcearrays(i_source,1,i,i)=',sourcearrays(i_source,1,xx,yy),'i=',xx,'j=',yy
             ! print *, '*********########*********sourcearrays(i_source,2,i,j)=',sourcearrays(i_source,2,xx,yy),'i=',xx,'j=',yy
             ! enddo
             ! enddo
             ! endif 
                
  do irec_local = 1,nrecloc
      ! if(it_tmp==7200-100) then
      !        print *, '*********########********* test myrank'
      !        print *, '*********########********* myrank=',myrank  
      ! endif
     i_source=irec_local
              
    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local) 
    if (ispec_is_elastic(ispec)) then
      ! add source array
      if (P_SV) then
        ! P-SV waves
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
             !lucas debug
             !if(mod(it_tmp,50)==0) then
             !print *, '********#######*****accel_elastic(2,iglob), before = ',accel_elastic(2,iglob), 'and it_tmp=',it_tmp
             !endif

          !  accel_elastic(1,iglob) = accel_elastic(1,iglob) + real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
          !                              source_adjoint(irec_local,it_tmp,1),kind=CUSTOM_REAL)
          !  accel_elastic(2,iglob) = accel_elastic(2,iglob) + real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
          !                              source_adjoint(irec_local,it_tmp,2),kind=CUSTOM_REAL)
          if(i==5 .and. j==5 ) then
          ! for plane wave
          sourcearrays_lucas(1,i,j)=sin(3.1415926-anglesource(i_source))
          sourcearrays_lucas(2,i,j)=-cos(3.1415926-anglesource(i_source))
          ! for point source
          !sourcearrays_lucas(1,i,j)=sin(anglesource(i_source))
          !sourcearrays_lucas(2,i,j)=-cos(anglesource(i_source))
          else
          sourcearrays_lucas(1,i,j)=0;
          sourcearrays_lucas(2,i,j)=0;
          endif
            ! lucas, the adjoint source, I use source_adjoint(irec_local,it_tmp,2) for vertical, and source_adjoint(irec_local,it_tmp,1) for horizontal

            !accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays_lucas(1,i,j)*source_adjoint(irec_local,it_tmp,1) ! lucas added angle here, sin(132.3*3.1415926/180.d0)
            !accel_elastic(2,iglob) = accel_elastic(2,iglob) + sourcearrays_lucas(2,i,j)*source_adjoint(irec_local,it_tmp,2) ! lucas added angle here, (-cos(132.3*3.1415926/180.d0))*

           !accel_elastic(1,iglob) = accel_elastic(1,iglob) + source_adjoint(irec_local,it_tmp,1)
                                   ! sin(3.1415926-anglesource(i_source)) ! lucas added angle here, sin(132.3*3.1415926/180.d0)
           !accel_elastic(2,iglob) = accel_elastic(2,iglob) + source_adjoint(irec_local,it_tmp,2)
                                   ! (-cos(3.1415926-anglesource(i_source)))  ! lucas (-cos(132.3*3.1415926/180.d0))

            
            accel_elastic(1,iglob) = accel_elastic(1,iglob) + real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)*&
                                     source_adjoint(irec_local,it_tmp,1) ,kind=CUSTOM_REAL) !*sin(100*3.1415926/180.d0) !lucas 155 for 25 incoming wave
            accel_elastic(2,iglob) = accel_elastic(2,iglob) + real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)*&
                                     source_adjoint(irec_local,it_tmp,2) ,kind=CUSTOM_REAL) !*(-cos(100*3.1415926/180.d0)) !lucas 155 for 25 incoming wave
           
              !lucas debug
              if(it_tmp==7200-100 .and. i==5 .and. j==5) then
              print *, '*********########********* in the loop'
              print *, '*********########*********sourcearrays_lucas(1,5,5)=',sourcearrays_lucas(1,5,5)  
              print *, '*********########*********sourcearrays_lucas(2,5,5)=',sourcearrays_lucas(2,5,5)
              !print *, '*********########*********xir_store_loc(irec_local,i)=',xir_store_loc(irec_local,i), 'i=', i
              !print *, '*********########*********gammar_store_loc(irec_local,j)=',gammar_store_loc(irec_local,j), 'j=', j
              print *, '*********########*********accel_elastic(1,iglob)=',accel_elastic(1,iglob)
              print *, '*********########*********accel_elastic(2,iglob)=',accel_elastic(2,iglob)
              print *, '*********########********* source_adjoint(irec_local,it_tmp,1)=',source_adjoint(irec_local,it_tmp,1)
              print *, '*********########********* source_adjoint(irec_local,it_tmp,2)=',source_adjoint(irec_local,it_tmp,2)
              endif


          enddo
        enddo

      else
        ! SH (membrane) wavescompute_forces_v
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            accel_elastic(1,iglob) = accel_elastic(1,iglob) +  real(xir_store_loc(irec_local,i)*&
                gammar_store_loc(irec_local,j)*source_adjoint(irec_local,it_tmp,1),kind=CUSTOM_REAL)

          enddo
        enddo
      endif
    endif ! if element is elastic
  enddo ! irec_local = 1,nrecloc

               !if(it_tmp==7200-100) then
               !print *, '*********########********* in the end'
               !print *, '*********########*********sourcearrays(i_source,1,5,5)=',sourcearrays(i_source,1,5,5)  
               !print *, '*********########*********sourcearrays(i_source,2,5,5)=',sourcearrays(i_source,2,5,5)
               !endif




  end subroutine compute_add_sources_viscoelastic_adjoint

 !-----------------lucas, CTD-SEM-------m2------------------------------------------------------------------------
 ! for viscoelastic solver for adjoint propagation wave field

  subroutine compute_add_sources_viscoelastic_adjoint_m2() 

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM !lucas added NDIM
  

  use specfem_par, only: P_SV,accel_elastic_m2, ispec_is_elastic,NSTEP,it, & !lucas, CTD-SEM
                         nrecloc,ispec_selected_rec_loc,ibool, & !myrank, & ! lucas added myrank 
                         source_adjoint,xir_store_loc,gammar_store_loc, anglesource ! sourcearrays
                         !lucas, xir_store_loc is the lagrange interpolant of receiver, see setup_source_receivers.f90, !lucas add anglesource and sourcearrays
                         
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec,i_source ! lucas added i_source=1
  integer :: it_tmp !,xx,yy
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearrays_lucas

  ! time step index
  it_tmp = NSTEP - it + 1
              i_source=1
             
             ! if(it_tmp==7200-100) then
             ! do yy=1,NGLLZ
             ! do xx=1,NGLLX
             ! print *, '*********########********* before the do'
             ! print *, '*********########********* myrank=',myrank
             ! print *, '*********########*********sourcearrays(i_source,1,i,i)=',sourcearrays(i_source,1,xx,yy),'i=',xx,'j=',yy
             ! print *, '*********########*********sourcearrays(i_source,2,i,j)=',sourcearrays(i_source,2,xx,yy),'i=',xx,'j=',yy
             ! enddo
             ! enddo
             ! endif 
                
  do irec_local = 1,nrecloc
       !if(it_tmp==7200-100) then
       !       print *, '*********########********* test myrank'
       !       print *, '*********########********* myrank=',myrank  
       !endif
     i_source=irec_local
              
    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local) 
    if (ispec_is_elastic(ispec)) then
      ! add source array
      if (P_SV) then
        ! P-SV waves
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

          if(i==5 .and. j==5 ) then
          ! for plane wave
          sourcearrays_lucas(1,i,j)=sin(3.1415926-anglesource(i_source))
          sourcearrays_lucas(2,i,j)=-cos(3.1415926-anglesource(i_source))
          ! for point source
          !sourcearrays_lucas(1,i,j)=sin(anglesource(i_source))
          !sourcearrays_lucas(2,i,j)=-cos(anglesource(i_source))
          else
          sourcearrays_lucas(1,i,j)=0;
          sourcearrays_lucas(2,i,j)=0;
          endif
            ! lucas, the adjoint source, I use source_adjoint(irec_local,it_tmp,2) for vertical, and source_adjoint(irec_local,it_tmp,1) for horizontal

            accel_elastic_m2(1,iglob) = accel_elastic_m2(1,iglob) + real(xir_store_loc(irec_local,i)*&
                gammar_store_loc(irec_local,j)*source_adjoint(irec_local,it_tmp,1),kind=CUSTOM_REAL) !*sin(100*3.1415926/180.d0)
            accel_elastic_m2(2,iglob) = accel_elastic_m2(2,iglob) + real(xir_store_loc(irec_local,i)*&
                gammar_store_loc(irec_local,j)*source_adjoint(irec_local,it_tmp,2),kind=CUSTOM_REAL) !*(-cos(100*3.1415926/180.d0))

          enddo
        enddo

      else
        ! SH (membrane) wavescompute_forces_v
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            
            accel_elastic_m2(1,iglob) = accel_elastic_m2(1,iglob) +  real(xir_store_loc(irec_local,i)*&
               gammar_store_loc(irec_local,j)*source_adjoint(irec_local,it_tmp,1),kind=CUSTOM_REAL)
            


          enddo
        enddo
      endif
    endif ! if element is elastic
  enddo ! irec_local = 1,nrecloc

               !if(it_tmp==7200-100) then
               !print *, '*********########********* in the end'
               !print *, '*********########*********sourcearrays(i_source,1,5,5)=',sourcearrays(i_source,1,5,5)  
               !print *, '*********########*********sourcearrays(i_source,2,5,5)=',sourcearrays(i_source,2,5,5)
               !endif




  end subroutine compute_add_sources_viscoelastic_adjoint_m2

 !-----------------lucas, CTD-SEM------m1-------------------------------------------------------------------------
 ! for viscoelastic solver for adjoint propagation wave field , for approximate Hessian

  subroutine compute_add_sources_viscoelastic_adjoint_m1() 

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM !lucas added NDIM
  

  use specfem_par, only: P_SV,accel_elastic_m1, ispec_is_elastic,NSTEP,it, & !lucas, CTD-SEM
                         nrecloc,ispec_selected_rec_loc,ibool, & !  myrank, & ! lucas added myrank
                         source_adjoint_m2,xir_store_loc,gammar_store_loc, anglesource !lucas, CTD-SEM
                         !lucas, xir_store_loc is the lagrange interpolant of receiver, see setup_source_receivers.f90, !lucas add anglesource and sourcearrays
                         
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec,i_source ! lucas added i_source=1
  integer :: it_tmp !,xx,yy
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearrays_lucas

  ! time step index
  it_tmp = NSTEP - it + 1
              i_source=1
             
             ! if(it_tmp==7200-100) then
             ! do yy=1,NGLLZ
             ! do xx=1,NGLLX
             ! print *, '*********########********* before the do'
             ! print *, '*********########********* myrank=',myrank
             ! print *, '*********########*********sourcearrays(i_source,1,i,i)=',sourcearrays(i_source,1,xx,yy),'i=',xx,'j=',yy
             ! print *, '*********########*********sourcearrays(i_source,2,i,j)=',sourcearrays(i_source,2,xx,yy),'i=',xx,'j=',yy
             ! enddo
             ! enddo
             ! endif 
                
  do irec_local = 1,nrecloc
       !if(it_tmp==7200-100) then
       !       print *, '*********########********* test myrank'
       !       print *, '*********########********* myrank=',myrank  
       !endif
     i_source=irec_local
              
    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local) 
    if (ispec_is_elastic(ispec)) then
      ! add source array
      if (P_SV) then
        ! P-SV waves
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

          if(i==5 .and. j==5 ) then
          ! for plane wave
          sourcearrays_lucas(1,i,j)=sin(3.1415926-anglesource(i_source))
          sourcearrays_lucas(2,i,j)=-cos(3.1415926-anglesource(i_source))
          ! for point source
          !sourcearrays_lucas(1,i,j)=sin(anglesource(i_source))
          !sourcearrays_lucas(2,i,j)=-cos(anglesource(i_source))
          else
          sourcearrays_lucas(1,i,j)=0;
          sourcearrays_lucas(2,i,j)=0;
          endif
            ! lucas, the adjoint source, I use source_adjoint_m2(irec_local,it_tmp,2) for vertical, and source_adjoint_m2(irec_local,it_tmp,1) for horizontal
            ! lucas, source_adjoint_m2 is the adjoint source computed by the waveforms of m2
            accel_elastic_m1(1,iglob) = accel_elastic_m1(1,iglob) + real(xir_store_loc(irec_local,i)*&
                gammar_store_loc(irec_local,j)*source_adjoint_m2(irec_local,it_tmp,1),kind=CUSTOM_REAL) !*sin(100*3.1415926/180.d0)
            accel_elastic_m1(2,iglob) = accel_elastic_m1(2,iglob) + real(xir_store_loc(irec_local,i)*&
                gammar_store_loc(irec_local,j)*source_adjoint_m2(irec_local,it_tmp,2),kind=CUSTOM_REAL) !*(-cos(100*3.1415926/180.d0))

          enddo
        enddo

      else
        ! SH (membrane) wavescompute_forces_v
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            
            accel_elastic_m1(1,iglob) = accel_elastic_m1(1,iglob) +  real(xir_store_loc(irec_local,i)*&
               gammar_store_loc(irec_local,j)*source_adjoint_m2(irec_local,it_tmp,1),kind=CUSTOM_REAL)
            


          enddo
        enddo
      endif
    endif ! if element is elastic
  enddo ! irec_local = 1,nrecloc

               !if(it_tmp==7200-100) then
               !print *, '*********########********* in the end'
               !print *, '*********########*********sourcearrays(i_source,1,5,5)=',sourcearrays(i_source,1,5,5)  
               !print *, '*********########*********sourcearrays(i_source,2,5,5)=',sourcearrays(i_source,2,5,5)
               !endif




  end subroutine compute_add_sources_viscoelastic_adjoint_m1


