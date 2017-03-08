module fault_module

    implicit none
    real(kind=8), parameter :: LAT2METER = 111133.84012073893 ! from clawpack.geoclaw.data
    integer :: nsubfaults, nevents
    real(kind=8), allocatable :: width(:), depth(:), slip(:), longitude(:), rupture_time(:), rise_time(:)
    real(kind=8), allocatable :: event_times(:)
    real(kind=8) :: center(2), theta, xcb(2)

contains

    subroutine load_fault(fname, nsubfaults)

        implicit none

        character*12, intent(in) :: fname
        integer, intent(in) :: nsubfaults

        integer :: i, j, k
        real(kind=8) :: input_line(12), xp1, yp1, xp2, yp2, total_width, swap

        allocate(width(nsubfaults))
        allocate(depth(nsubfaults))
        allocate(slip(nsubfaults))
        allocate(longitude(nsubfaults))
        allocate(rupture_time(nsubfaults))
        allocate(rise_time(nsubfaults))
        allocate(event_times(2*nsubfaults))

        call opendatafile(7, fname)

        read(7,*)
        read(7,*) !nsubfaults
        read(7,*)

        ! Read in subfaults
        xp1 = 1.d10
        xp2 = -1.d10
        write(*,*) 'Here 1'
        do i=1,nsubfaults
          read(7,*) input_line
          write(*,*) 'done reading'
          theta = input_line(2)/180.0*4.d0*datan(1.d0)
          width(i) = input_line(3)
          depth(i) = input_line(4)
          slip(i) = input_line(5)
          longitude(i) = input_line(9)
          rupture_time(i) = input_line(11)
          rise_time(i) = input_line(12)

          event_times(2*i-1) = rupture_time(i)
          event_times(2*i) = rupture_time(i) + rise_time(i)

          swap = longitude(i)*LAT2METER
          if (swap < xp1) then
            xp1 = swap
            yp1 = -depth(i)
          end if

          swap = longitude(i)*LAT2METER + dcos(theta)*width(i)
          if (xp2 < swap) then
            xp2 = swap
            yp2 = -depth(i) - dsin(theta)*width(i)
          end if
        end do
        write(*,*) 'Here 2'


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
        write(*,*) 'Here 3'

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
        write(*,*) 'Here 4'

        ! Compute spatial fault information
        center(1) = 0.5d0*(xp1 + xp2)
        center(2) = 0.5d0*(yp1 + yp2)
        total_width = dsqrt((xp2-xp1)**2 + (yp2-yp1)**2)

        xcb(1) = center(1) - 0.5d0*total_width
        xcb(2) = center(1) + 0.5d0*total_width

        close(7)
        write(*,*) 'Here 5'

    end subroutine load_fault

end module fault_module
