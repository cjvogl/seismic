!   ==================
    subroutine setprob
!   ==================

    implicit none

    character*12 fname
    integer iunit

    REAL (kind=8) :: center(2), theta, xcb(2), mindepth
    common /fault/  center, theta, xcb, mindepth

    real (kind=8) :: scaling
    common /water/  scaling

    real (kind=8) :: width
!
!
      iunit = 7
      fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
      call opendatafile(iunit, fname)


!
      read(7,*) width ! domain_depth
      read(7,*) width ! domain width
      read(7,*) center(1)
      read(7,*) width
      read(7,*) theta
      read(7,*) center(2)
      read(7,*) scaling ! water depth
      read(7,*) scaling

      center(2) = -center(2)
      xcb(1) = center(1) - 0.5*width
      xcb(2) = center(1) + 0.5*width

      mindepth = dmin1(dabs(center(2) - 0.5*width*dsin(theta)), &
                      dabs(center(2) + 0.5*width*dsin(theta)))

    return
    end
