subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

!   set slip before call to step2

    use fault_module, only: center, xcb, nsubfaults, subfaults
    use fault_module, only: nevents, event_times

    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    real(kind=8) :: fault_zshift
    common /mapping/ fault_zshift

    integer :: i, j, k
    real(kind=8) :: xcell, ycell, xpcell, ypcell

    aux(13,:,:) = 0.d0
    if (t <= event_times(nevents)) then

      do j=1-mbc,my+mbc
        ycell = ylower + (j-0.5d0)*dy
        if (abs(ycell - 0.5d0*dy + fault_zshift - center(2)) < 0.5d0*dy) then

          do i=1-mbc,mx+mbc
            xcell = xlower + (i-0.5d0)*dx
            if (xcb(1)-1.d-10 <= xcell - 0.5d0*dx .and. xcell + 0.5d0*dx <= xcb(2)+1.d-10) then
              ! find which subfault this cell center lies in and apply slip
              do k=1,nsubfaults
                if (subfaults(k)%xcb(1) <= xcell .and. &
                    xcell <= subfaults(k)%xcb(2) .and. &
                    subfaults(k)%rupture_time <= t .and. &
                    t <= subfaults(k)%rupture_time + subfaults(k)%rise_time) then
                  aux(13,i,j) = subfaults(k)%slip/subfaults(k)%rise_time
                  exit
                end if
              end do

            end if
          end do

        end if
      end do

    end if

end subroutine b4step2
