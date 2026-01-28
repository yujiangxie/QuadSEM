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

  subroutine create_color_image(i_field,plot_b_wavefield_only) !lucas, plot_b_wavefield_only=.false.

! display a given field as a red and blue color JPEG image

! to display the snapshots : display image*.jpg

  use constants, only: TINYVAL,VERYTINYVAL,HUGEVAL,STABILITY_THRESHOLD,OUTPUT_FILES,IMAIN

  use specfem_par, only: myrank,it,NSOURCES,P_SV,nrec,CTD_SEM,Full_Hessian_by_Wavefield_Stored,compute_appro_Hessian !lucas, CTD-SEM

  use specfem_par_movie, only: image_color_data,iglob_image_color,NX_IMAGE_color,NZ_IMAGE_color, &
    isnapshot_number,cutsnaps,image_color_vp_display,image_color_vp_display_m2,image_color_vp_display_m1, & !lucas, CTD-SEM
    image_color_data_m2,image_color_data_m1, & !lucas, CTD-SEM
    USE_SNAPSHOT_NUMBER_IN_FILENAME,POWER_DISPLAY_COLOR, &
    DRAW_SOURCES_AND_RECEIVERS, &
    ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver, &
    USE_CONSTANT_MAX_AMPLITUDE,CONSTANT_MAX_AMPLITUDE_TO_USE,SIMULATION_TYPE

  implicit none

  integer :: i_field
  logical :: plot_b_wavefield_only

  ! local parameters
  integer :: i

  ! for the JPEG library
  character(len=1), dimension(3,NX_IMAGE_color,NZ_IMAGE_color) :: JPEG_raw_image,JPEG_raw_image_m2,JPEG_raw_image_m1 !lucas, CTD-SEM
  integer :: ix,iy,R,G,B
  double precision :: amplitude_max,normalized_value,vpmin,vpmax,x1,amplitude_max_m2,normalized_value_m2,vpmin_m2,vpmax_m2 !lucas, CTD-SEM
  double precision :: amplitude_max_m1,vpmin_m1,vpmax_m1,normalized_value_m1 !lucas, CTD-SEM
  character(len=150) :: filename, filename_m2, filename_m11 !lucas, CTD-SEM
  logical :: do_warning, do_warning_m2,do_warning_m1 !lucas, CTD-SEM

  ! size of cross and square in pixels drawn to represent the source and the receivers in JPEG pictures
  integer :: half_width_cross, thickness_cross, half_size_square

! make the size of the source and receiver symbols depend on the size of the picture
! using a rule of thumb
  thickness_cross = 1
  if (NX_IMAGE_color > 2000 .or. NZ_IMAGE_color > 2000) then
    half_width_cross = 6
    half_size_square = 4
  else if (NX_IMAGE_color <= 100 .or. NZ_IMAGE_color <= 100) then
    half_width_cross = 2
    half_size_square = 1
  else if (NX_IMAGE_color <= 250 .or. NZ_IMAGE_color <= 250) then
    half_width_cross = 3
    half_size_square = 2
  else
    half_width_cross = 5
    half_size_square = 3
  endif

! open the image file
! slightly change the beginning of the file name depending if we use the time step of the image number, to avoid confusion
  if (USE_SNAPSHOT_NUMBER_IN_FILENAME) then !lucas, USE_SNAPSHOT_NUMBER_IN_FILENAME=.false.
    isnapshot_number = isnapshot_number + 1
    if (i_field == 1 .and. SIMULATION_TYPE == 1) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'forward_img',isnapshot_number,'.jpg'
    else if (i_field == 1 .and. SIMULATION_TYPE == 3 .and. .not. plot_b_wavefield_only) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'adjoint_img',isnapshot_number,'.jpg'
    else
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'b_forward_img',isnapshot_number,'.jpg'
    endif
  else
    if (i_field == 1 .and. SIMULATION_TYPE == 1) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'forward_image',it,'.jpg'

      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas1, CTD-SEM
      write(filename_m2,"(a,i7.7,a)") trim(OUTPUT_FILES)//'forward_image_m2_',it,'.jpg'
      endif

    else if (i_field == 1 .and. SIMULATION_TYPE == 3 .and. .not. plot_b_wavefield_only) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'adjoint_image',it,'.jpg'
      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas1.1, CTD-SEM
       write(filename_m2,"(a,i7.7,a)") trim(OUTPUT_FILES)//'adjoint_image_m2_',it,'.jpg'
       if(compute_appro_Hessian) then
          write(filename_m11,"(a,i7.7,a)") trim(OUTPUT_FILES)//'adjoint_image_m1_adjsrc-m2_',it,'.jpg' !lucas, approximate Hessian only need one here
       endif
      endif
    else
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'b_forward_image',it,'.jpg'
      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas1.2, CTD-SEM
      write(filename_m2,"(a,i7.7,a)") trim(OUTPUT_FILES)//'b_forward_image_m2_',it,'.jpg'
      endif
    endif
  endif

! compute maximum amplitude
  if (.not. USE_CONSTANT_MAX_AMPLITUDE) then !lucas, USE_CONSTANT_MAX_AMPLITUDE=.false.
    amplitude_max = maxval(abs(image_color_data))
    if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas2, CTD-SEM
     amplitude_max_m2 = maxval(abs(image_color_data_m2))
     if(compute_appro_Hessian) then
       amplitude_max_m1 = maxval(abs(image_color_data_m1))
     endif
    endif 

    !amplitude_max =2.0402854961076855d-9 !lucas----------m+dm-------------------------------------------------for one 2D model only (my Hessian paper)
    !amplitude_max = max(amplitude_max,1.2277411309824515d-9) !lucas-----------------------------------------------------------for one 2D model only (my Hessian paper)
    print *,'min and max amplitude (Lucas) for image_color_data(i,j)= ',minval(image_color_data), amplitude_max
    if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas3, CTD-SEM
     print *,'min and max amplitude (Lucas) for image_color_data_m2(i,j)= ',minval(image_color_data_m2), amplitude_max_m2
     if(compute_appro_Hessian) then
       print *,'min and max amplitude (Lucas) for image_color_data_m1(i,j)= ',minval(image_color_data_m1), amplitude_max_m1
     endif
    endif

  else 
    amplitude_max = CONSTANT_MAX_AMPLITUDE_TO_USE
!   in case of a pre-defined and constant maximum, truncate all values that are outside that constant range
    where(image_color_data > +CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = +CONSTANT_MAX_AMPLITUDE_TO_USE
    where(image_color_data < -CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = -CONSTANT_MAX_AMPLITUDE_TO_USE
    if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas4, CTD-SEM
     amplitude_max_m2 = CONSTANT_MAX_AMPLITUDE_TO_USE
     where(image_color_data_m2 > +CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data_m2 = +CONSTANT_MAX_AMPLITUDE_TO_USE
     where(image_color_data_m2 < -CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data_m2 = -CONSTANT_MAX_AMPLITUDE_TO_USE
     if(compute_appro_Hessian) then
        amplitude_max_m1 = CONSTANT_MAX_AMPLITUDE_TO_USE
        where(image_color_data_m1 > +CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data_m1 = +CONSTANT_MAX_AMPLITUDE_TO_USE
        where(image_color_data_m1 < -CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data_m1 = -CONSTANT_MAX_AMPLITUDE_TO_USE
     endif
    endif

  endif
  ! user output
  write(IMAIN,*) 'Color image maximum amplitude = ',amplitude_max
  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas5, CTD-SEM
   write(IMAIN,*) 'Color image maximum amplitude for m2 in CTD-SEM = ',amplitude_max_m2
   if(compute_appro_Hessian) then
     write(IMAIN,*) 'Color image maximum amplitude for m1 in CTD-SEM = ',amplitude_max_m1
   endif
  endif

! this trick checks for NaN (Not a Number), which is not even equal to itself
  if (amplitude_max > STABILITY_THRESHOLD .or. amplitude_max < 0 .or. amplitude_max /= amplitude_max) then
    print *,'Warning: failed creating color image, maximum value of amplitude in image color is invalid'
    print *,'amplitude max = ',amplitude_max,' with threshold at ', STABILITY_THRESHOLD
    print *,'Please check your simulation setup...'
    call exit_MPI(myrank,'code became unstable and blew up (image_color_data)')
  endif

  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas6, CTD-SEM
   if (amplitude_max_m2 > STABILITY_THRESHOLD .or. amplitude_max_m2 < 0 .or. amplitude_max_m2 /= amplitude_max_m2) then
    print *,'Warning: failed creating color image, maximum value of amplitude in image color is invalid'
    print *,'amplitude max_m2 = ',amplitude_max_m2,' with threshold at ', STABILITY_THRESHOLD
    print *,'Please check your simulation setup...'
    call exit_MPI(myrank,'code became unstable and blew up (image_color_data)')
   endif
   if(compute_appro_Hessian) then !----
    if (amplitude_max_m1 > STABILITY_THRESHOLD .or. amplitude_max_m1 < 0 .or. amplitude_max_m1 /= amplitude_max_m1) then
     print *,'Warning: failed creating color image, maximum value of amplitude in image color is invalid'
     print *,'amplitude max_m1 = ',amplitude_max_m1,' with threshold at ', STABILITY_THRESHOLD
     print *,'Please check your simulation setup...'
     call exit_MPI(myrank,'code became unstable and blew up (image_color_data)')
    endif
   endif !---
  endif


  do_warning = .false.
  do_warning_m2 = .false.
  do_warning_m1 = .false.

  vpmin = HUGEVAL
  vpmax = TINYVAL
  
  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas6_1, CTD-SEM
   vpmin_m2 = HUGEVAL
   vpmax_m2 = TINYVAL
   if(compute_appro_Hessian) then
    vpmin_m1 = HUGEVAL
    vpmax_m1 = TINYVAL
   endif
  endif
  
  do iy= 1,NZ_IMAGE_color
    do ix= 1,NX_IMAGE_color
! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue later
      if (iglob_image_color(ix,iy) > -1 .and. image_color_vp_display(ix,iy) >= 0) then
        vpmin = min(vpmin,image_color_vp_display(ix,iy)) !lucas image_color_vp_display=vp values.
        vpmax = max(vpmax,image_color_vp_display(ix,iy))
      endif

      if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas6_2, CTD-SEM
       if (iglob_image_color(ix,iy) > -1 .and. image_color_vp_display_m2(ix,iy) >= 0) then
        vpmin_m2 = min(vpmin_m2,image_color_vp_display_m2(ix,iy)) !lucas image_color_vp_display=vp values.
        vpmax_m2 = max(vpmax_m2,image_color_vp_display_m2(ix,iy))
       endif
       if(compute_appro_Hessian) then !---
        if (iglob_image_color(ix,iy) > -1 .and. image_color_vp_display_m1(ix,iy) >= 0) then
          vpmin_m1 = min(vpmin_m1,image_color_vp_display_m1(ix,iy)) !lucas image_color_vp_display=vp values.
          vpmax_m1 = max(vpmax_m1,image_color_vp_display_m1(ix,iy))
        endif
       endif !----
      endif

    enddo
  enddo

! in the image format, the image starts in the upper-left corner
  do iy=NZ_IMAGE_color,1,-1
    do ix= 1,NX_IMAGE_color

! check if pixel is defined or not (can be above topography for instance)
      if (iglob_image_color(ix,iy) == -1) then

! use white to display undefined region above topography to avoid visual confusion with a water layer
        R = 255
        G = 255
        B = 255

! suppress small amplitudes considered as noise and display the background velocity model instead
      else if (abs(image_color_data(ix,iy)) < amplitude_max * cutsnaps) then

! use P velocity model as background where amplitude is negligible
        if ((P_SV) .and. ((vpmax-vpmin)/max(vpmin, TINYVAL) > 0.02d0)) then
          x1 = (image_color_vp_display(ix,iy)-vpmin)/(vpmax-vpmin) 
        else
          x1 = 0.5d0 
        endif
          
! rescale to avoid very dark gray levels
        x1 = x1*0.7 + 0.2
        if (x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin, dark gray = vpmax
        x1 = 1.d0 - x1

! map to [0,255]
        x1 = x1 * 255.d0

        R = nint(x1)
        if (R < 0) R = 0
        if (R > 255) R = 255
        G = R 
        B = R

! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue
        if (image_color_vp_display(ix,iy) < 0) then
! use light blue to display water
!!!!!!    R = 204
!!!!!!    G = 255
!!!!!!    B = 255
          R = 135 !!! LightSkyBlue
          G = 206
          B = 250
        endif

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        if (amplitude_max >= VERYTINYVAL) then
           normalized_value = image_color_data(ix,iy) / amplitude_max
          !normalized_value = image_color_data(ix,iy) ! lucas, no norm----------------------- not work yet
        else
          normalized_value = image_color_data(ix,iy) / VERYTINYVAL
        endif
        ! check value (isNaN)
        if (normalized_value /= normalized_value) then
          ! will be set to zero
          normalized_value = 0.d0
          do_warning = .true.
        endif

! suppress values outside of [-1:+1]
        if (normalized_value < -1.d0) normalized_value = -1.d0 
        if (normalized_value > 1.d0) normalized_value = 1.d0 

! use red if positive value, blue if negative, no green  ! lucas, wavefiled color
        if (normalized_value >= 0.d0) then
          R = nint(255.d0*normalized_value**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY_COLOR)
        endif

      endif

! for JPEG
     JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
     JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
     JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)

     !lucas7 -------------------------- CTD-SEM-----start---------------------------------------------------------
     if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then

     ! check if pixel is defined or not (can be above topography for instance)
      if (iglob_image_color(ix,iy) == -1) then

! use white to display undefined region above topography to avoid visual confusion with a water layer
        R = 255
        G = 255
        B = 255

! suppress small amplitudes considered as noise and display the background velocity model instead
      else if (abs(image_color_data_m2(ix,iy)) < amplitude_max_m2 * cutsnaps) then

! use P velocity model as background where amplitude is negligible
        if ((P_SV) .and. ((vpmax_m2-vpmin_m2)/max(vpmin_m2, TINYVAL) > 0.02d0)) then
          x1 = (image_color_vp_display_m2(ix,iy)-vpmin_m2)/(vpmax_m2-vpmin_m2) 
        else
          x1 = 0.5d0 
        endif
          
! rescale to avoid very dark gray levels
        x1 = x1*0.7 + 0.2
        if (x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin_m2, dark gray = vpmax_m2
        x1 = 1.d0 - x1

! map to [0,255]
        x1 = x1 * 255.d0

        R = nint(x1)
        if (R < 0) R = 0
        if (R > 255) R = 255
        G = R 
        B = R

! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue
        if (image_color_vp_display_m2(ix,iy) < 0) then
! use light blue to display water
!!!!!!    R = 204
!!!!!!    G = 255
!!!!!!    B = 255
          R = 135 !!! LightSkyBlue
          G = 206
          B = 250
        endif

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        if (amplitude_max_m2 >= VERYTINYVAL) then
           normalized_value_m2 = image_color_data_m2(ix,iy) / amplitude_max_m2
        else
          normalized_value_m2 = image_color_data_m2(ix,iy) / VERYTINYVAL
        endif
        ! check value (isNaN)
        if (normalized_value_m2 /= normalized_value_m2) then
          ! will be set to zero
          normalized_value_m2 = 0.d0
          do_warning_m2 = .true.
        endif

! suppress values outside of [-1:+1]
        if (normalized_value_m2 < -1.d0) normalized_value_m2 = -1.d0 
        if (normalized_value_m2 > 1.d0) normalized_value_m2 = 1.d0 

! use red if positive value, blue if negative, no green  ! lucas, wavefiled color
        if (normalized_value_m2 >= 0.d0) then
          R = nint(255.d0*normalized_value_m2**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(normalized_value_m2)**POWER_DISPLAY_COLOR)
        endif

      endif

! for JPEG
     JPEG_raw_image_m2(1,ix,NZ_IMAGE_color-iy+1) = char(R)
     JPEG_raw_image_m2(2,ix,NZ_IMAGE_color-iy+1) = char(G)
     JPEG_raw_image_m2(3,ix,NZ_IMAGE_color-iy+1) = char(B)
      
     ! when under the if condition of CTD_SEM or the stored method
     if(compute_appro_Hessian) then !#######################################################################################
             ! check if pixel is defined or not (can be above topography for instance)
      if (iglob_image_color(ix,iy) == -1) then

! use white to display undefined region above topography to avoid visual confusion with a water layer
        R = 255
        G = 255
        B = 255

! suppress small amplitudes considered as noise and display the background velocity model instead
      else if (abs(image_color_data_m1(ix,iy)) < amplitude_max_m1 * cutsnaps) then

! use P velocity model as background where amplitude is negligible
        if ((P_SV) .and. ((vpmax_m1-vpmin_m1)/max(vpmin_m1, TINYVAL) > 0.02d0)) then
          x1 = (image_color_vp_display_m1(ix,iy)-vpmin_m1)/(vpmax_m1-vpmin_m1) 
        else
          x1 = 0.5d0 
        endif
          
! rescale to avoid very dark gray levels
        x1 = x1*0.7 + 0.2
        if (x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin_m1, dark gray = vpmax_m1
        x1 = 1.d0 - x1

! map to [0,255]
        x1 = x1 * 255.d0

        R = nint(x1)
        if (R < 0) R = 0
        if (R > 255) R = 255
        G = R 
        B = R

! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue
        if (image_color_vp_display_m1(ix,iy) < 0) then
! use light blue to display water
!!!!!!    R = 204
!!!!!!    G = 255
!!!!!!    B = 255
          R = 135 !!! LightSkyBlue
          G = 206
          B = 250
        endif

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        if (amplitude_max_m1 >= VERYTINYVAL) then
           normalized_value_m1 = image_color_data_m1(ix,iy) / amplitude_max_m1
        else
          normalized_value_m1 = image_color_data_m1(ix,iy) / VERYTINYVAL
        endif
        ! check value (isNaN)
        if (normalized_value_m1 /= normalized_value_m1) then
          ! will be set to zero
          normalized_value_m1 = 0.d0
          do_warning_m1 = .true.
        endif

! suppress values outside of [-1:+1]
        if (normalized_value_m1 < -1.d0) normalized_value_m1 = -1.d0 
        if (normalized_value_m1 > 1.d0) normalized_value_m1 = 1.d0 

! use red if positive value, blue if negative, no green  ! lucas, wavefiled color
        if (normalized_value_m1 >= 0.d0) then
          R = nint(255.d0*normalized_value_m1**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(normalized_value_m1)**POWER_DISPLAY_COLOR)
        endif

      endif

! for JPEG
     JPEG_raw_image_m1(1,ix,NZ_IMAGE_color-iy+1) = char(R)
     JPEG_raw_image_m1(2,ix,NZ_IMAGE_color-iy+1) = char(G)
     JPEG_raw_image_m1(3,ix,NZ_IMAGE_color-iy+1) = char(B)
     endif !######################################################################################################
     endif !end lucas7, CTD-SEM 

     !lucas7-------------------------------CTD-SEM---end---------------------------------------------------------------
    enddo
  enddo
  

  if (do_warning) then
    print *,'Warning: normalized value is invalid, likely due to unstable simulation! process rank is ',myrank
    print *,'Warning: normalized_value = ',normalized_value,' will be set to zero'
  endif

  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas8, CTD-SEM
   if (do_warning_m2) then
    print *,'Warning: normalized value_m2 is invalid, likely due to unstable simulation! process rank is ',myrank
    print *,'Warning: normalized_value_m2 = ',normalized_value_m2,' will be set to zero'
   endif
   if(compute_appro_Hessian) then!---
    if (do_warning_m1) then
     print *,'Warning: normalized value_m1 is invalid, likely due to unstable simulation! process rank is ',myrank
     print *,'Warning: normalized_value_m1 = ',normalized_value_m1,' will be set to zero'
    endif
   endif !---
  endif

!
!----  draw position of the sources and receivers
!
  if (DRAW_SOURCES_AND_RECEIVERS) then  !lucas, DRAW_SOURCES_AND_RECEIVERS=.true.

! draw position of the sources with orange crosses
    do i = 1,NSOURCES

! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_source(i) - half_width_cross,1), min(iy_image_color_source(i) + half_width_cross,NZ_IMAGE_color)
        do ix = max(ix_image_color_source(i) - thickness_cross,1), min(ix_image_color_source(i) + thickness_cross,NX_IMAGE_color)
! use orange color
          R = 255
          G = 157
          B = 0
! for JPEG
          JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)
          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas9, CTD-SEM
           JPEG_raw_image_m2(1,ix,NZ_IMAGE_color-iy+1) = char(R)
           JPEG_raw_image_m2(2,ix,NZ_IMAGE_color-iy+1) = char(G)
           JPEG_raw_image_m2(3,ix,NZ_IMAGE_color-iy+1) = char(B)
           if(compute_appro_Hessian) then !---
             JPEG_raw_image_m1(1,ix,NZ_IMAGE_color-iy+1) = char(R)
             JPEG_raw_image_m1(2,ix,NZ_IMAGE_color-iy+1) = char(G)
             JPEG_raw_image_m1(3,ix,NZ_IMAGE_color-iy+1) = char(B)
           endif !---
          endif

        enddo
      enddo

! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_source(i) - thickness_cross,1), min(iy_image_color_source(i) + thickness_cross,NZ_IMAGE_color)
        do ix = max(ix_image_color_source(i) - half_width_cross,1), min(ix_image_color_source(i) + half_width_cross,NX_IMAGE_color)
! use orange color
          R = 255
          G = 157
          B = 0
! for JPEG
          JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)
          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas10, CTD-SEM
           JPEG_raw_image_m2(1,ix,NZ_IMAGE_color-iy+1) = char(R)
           JPEG_raw_image_m2(2,ix,NZ_IMAGE_color-iy+1) = char(G)
           JPEG_raw_image_m2(3,ix,NZ_IMAGE_color-iy+1) = char(B)
           if(compute_appro_Hessian) then !---
             JPEG_raw_image_m1(1,ix,NZ_IMAGE_color-iy+1) = char(R)
             JPEG_raw_image_m1(2,ix,NZ_IMAGE_color-iy+1) = char(G)
             JPEG_raw_image_m1(3,ix,NZ_IMAGE_color-iy+1) = char(B)
           endif !---
          endif 

        enddo
      enddo

    enddo

! draw position of the receivers with green squares
    do i = 1,nrec
! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_receiver(i) - half_size_square,1), &
                                          min(iy_image_color_receiver(i) + half_size_square,NZ_IMAGE_color)
        do ix = max(ix_image_color_receiver(i) - half_size_square,1), &
                                          min(ix_image_color_receiver(i) + half_size_square,NX_IMAGE_color)
! use dark green color
          R = 30
          G = 180
          B = 60
! for JPEG
          JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)
          if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas11, CTD-SEM
           JPEG_raw_image_m2(1,ix,NZ_IMAGE_color-iy+1) = char(R)
           JPEG_raw_image_m2(2,ix,NZ_IMAGE_color-iy+1) = char(G)
           JPEG_raw_image_m2(3,ix,NZ_IMAGE_color-iy+1) = char(B)
           if(compute_appro_Hessian) then !---
             JPEG_raw_image_m1(1,ix,NZ_IMAGE_color-iy+1) = char(R)
             JPEG_raw_image_m1(2,ix,NZ_IMAGE_color-iy+1) = char(G)
             JPEG_raw_image_m1(3,ix,NZ_IMAGE_color-iy+1) = char(B)
           endif !---
          endif


        enddo
      enddo
    enddo

  endif

! for JPEG
  call write_jpeg_image(JPEG_raw_image,NX_IMAGE_color,NZ_IMAGE_color,filename)
  if(CTD_SEM .or. Full_Hessian_by_Wavefield_Stored) then !lucas12, CTD-SEM
    call write_jpeg_image(JPEG_raw_image_m2,NX_IMAGE_color,NZ_IMAGE_color,filename_m2)
    if(compute_appro_Hessian) then !---
      call write_jpeg_image(JPEG_raw_image_m1,NX_IMAGE_color,NZ_IMAGE_color,filename_m11)
    endif !---
  endif 


  end subroutine create_color_image

