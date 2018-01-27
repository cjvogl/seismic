!   ==================
    subroutine setprob
!   ==================

    use fault_module, only: load_fault, center
    use amr_module, only: xlower,xupper,ylower

    implicit none

    character*12 fname
    integer iunit

    real(kind=8) :: rho1,amu1,alam1,rho2,amu2,alam2
    real(kind=8) :: ABLdepth, ABLxpos(2), ABLypos

    common /comaux/ rho1,amu1,alam1,rho2,amu2,alam2
    common /ablparam/ ABLdepth, ABLxpos, ABLypos

!
!     # fault parameters:
!
      fname = 'fault.data'
      call load_fault(fname)

!     # Set the material parameters for the elasticity equations
!
!
      iunit = 7
      fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
      call opendatafile(iunit, fname)

!
!     # Piecewise constant medium 
!     # Material parameters

      read(7,*) rho1
      read(7,*) alam1
      read(7,*) amu1

      read(7,*) rho2
      read(7,*) alam2
      read(7,*) amu2

      read(7,*) ABLdepth

!     # Compute ABL position
      ABLxpos(1) = xlower + ABLdepth
      ABLxpos(2) = xupper - ABLdepth
      ABLypos = ylower + ABLdepth
      write(6,*) '+++ ABLdepth, ABLypos: ',ABLdepth, ABLypos
      write(6,*) '+++ ABLxpos: ',ABLxpos

      return

    end subroutine setprob
