module fault_module

    implicit none
    real(kind=8), parameter :: LAT2METER = 111133.84012073893 ! from clawpack.geoclaw.data
    integer :: nsubfaults, nevents
    type subfault
      real(kind=8) :: width, depth, slip, longitude, rupture_time, rise_time
      real(kind=8) :: xcb(2), ycb(2)
    end type subfault
    type(subfault), allocatable :: subfaults(:)
    real(kind=8), allocatable :: event_times(:)
    real(kind=8) :: center(2), theta, xcb(2), ycb(2)

contains

    subroutine load_fault(fname)

        implicit none

        character*12, intent(in) :: fname

        integer :: i, j, k
        real(kind=8) :: input_line(12), xp1, yp1, xp2, yp2, total_width, swap

        call opendatafile(7, fname)

        read(7,*)
        read(7,*) nsubfaults
        read(7,*)

        allocate(subfaults(nsubfaults))
        allocate(event_times(2*nsubfaults))

        ! Read in subfaults
        xp1 = 1.d10
        xp2 = -1.d10
        do i=1,nsubfaults
          read(7,*) input_line
          theta = input_line(2)/180.0*4.d0*datan(1.d0)
          subfaults(i)%width = input_line(3)
          subfaults(i)%depth = input_line(4)
          subfaults(i)%slip = input_line(5)
          subfaults(i)%longitude = input_line(9)
          subfaults(i)%rupture_time = input_line(11)
          subfaults(i)%rise_time = input_line(12)

          event_times(2*i-1) = subfaults(i)%rupture_time
          event_times(2*i) = subfaults(i)%rupture_time + subfaults(i)%rise_time

          swap = subfaults(i)%longitude*LAT2METER
          if (swap < xp1) then
            xp1 = swap
            yp1 = -subfaults(i)%depth
          end if

          swap = subfaults(i)%longitude*LAT2METER + dcos(theta)*subfaults(i)%width
          if (xp2 < swap) then
            xp2 = swap
            yp2 = -subfaults(i)%depth - dsin(theta)*subfaults(i)%width
          end if
        end do
        close(7)

        ! Remove trivial event times
        nevents = 2*nsubfaults
        i = 1
        do while (i <= nevents)
          if (dabs(event_times(i)) < 1.d-10) then
            do j=i,nevents-1
              event_times(j) = event_times(j+1)
            end do
            event_times(nevents) = 0.d0
            nevents = nevents - 1
          else
            i = i+1
          end if
        end do

        ! Remove duplicate event times
        do while(i <= nevents-1)
          j = i+1
          do while(j <= nevents)
            if (dabs(event_times(i)-event_times(j)) < 1.d-10) then
              do k=j,nevents-1
                event_times(k) = event_times(k+1)
              end do
              event_times(nevents) = 0.d0
              nevents = nevents - 1
            else
              j = j+1
            end if

          end do
          i = i+1
        end do

        ! Sort event times
        do i=1,nevents-1
          do j=1,nevents-i
            if (event_times(j) > event_times(j+1)) then
              swap = event_times(j+1)
              event_times(j+1) = event_times(j)
              event_times(j) = swap
            end if
          end do
        end do

        ! Compute spatial fault information
        center(1) = 0.5d0*(xp1 + xp2)
        center(2) = 0.5d0*(yp1 + yp2)
        total_width = dsqrt((xp2-xp1)**2 + (yp2-yp1)**2)

        xcb(1) = center(1)
        xcb(2) = center(1)
        ycb(1) = center(2) - 0.5d0*total_width
        ycb(2) = center(2) + 0.5d0*total_width
        write(6,*) 'xcb: ',xcb(1),xcb(2)
        write(6,*) 'ycb: ',ycb(1),ycb(2)

        ! Determine subfault location in computational domain
        print *, "*******WARNING*******"
        print *, "Subfaults are assumed to be specified in top-center coords"
        print *, "*********************"
        do i=1,nsubfaults
          ! should be equal for vertical fault:
          subfaults(i)%xcb(1) = xcb(1) 
          subfaults(i)%xcb(2) = xcb(2)
          ! for vertical fault:
          subfaults(i)%ycb(2) = ycb(2) - (i-1)*subfaults(i)%width
          subfaults(i)%ycb(1) = ycb(2) - i*subfaults(i)%width
          write(6,*) '+++ subfault i, ycb1, ybc2: ', &
                i,subfaults(i)%ycb(1),subfaults(i)%ycb(2)
        end do

    end subroutine load_fault

end module fault_module
