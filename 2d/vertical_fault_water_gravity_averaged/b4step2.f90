subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

!   set slip before call to step2

    use fault_module, only: center, ycb, nsubfaults, subfaults
    use fault_module, only: nevents, event_times

    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    !real(kind=8) :: fault_zshift
    !common /mapping/ fault_zshift

    integer :: i, j, k
    real(kind=8) :: xcell, ycell, xpcell, ypcell, slip_per_second

    aux(6,:,:) = 0.d0   ! slip component for non-mapped Cartesian grid
    if (t <= event_times(nevents)) then

    !fault_zshift = 0.d0

    ! for a vertical fault (dip=90):

      do i=1-mbc,mx+mbc
        xcell = xlower + (i-0.5d0)*dx
        if (abs(xcell - 0.5d0*dx - center(1)) < 0.5d0*dx) then

          do j=1-mbc,my+mbc
            ycell = ylower + (j-0.5d0)*dy
            if (ycb(1)-1.d-10 <= ycell - 0.5d0*dy .and. &
                ycell + 0.5d0*dy <= ycb(2)+1.d-10) then
              ! find which subfault this cell center lies in and apply slip
!             write(6,*) '+++ along the fault, t, ycell = ',t, ycell
              do k=1,nsubfaults
!                     write(6,602) k, subfaults(k)%ycb(1), subfaults(k)%ycb(2)
 602                  format('+++ k, ycb1, ycb2: ', i3, 2f10.2)
                if (subfaults(k)%ycb(1) <= ycell .and. &
                    ycell <= subfaults(k)%ycb(2) .and. &
                    subfaults(k)%rupture_time <= t .and. &
                    t <= subfaults(k)%rupture_time &
                         + 2.d0*subfaults(k)%rise_time) then

                      slip_per_second = subfaults(k)%slip &
                            / (2.d0*subfaults(k)%rise_time)
                      aux(6,i,j) = 2.d0*slip_per_second ! jump in v

!                     write(6,601) xcell,ycell, dt, aux(6,i,j)
 601                  format('+++ xcell, ycell, dt, slip: ', 2f10.2,2f8.4)

                      exit
                    end if
              end do

            end if
          end do

        end if
      end do

    end if

end subroutine b4step2
