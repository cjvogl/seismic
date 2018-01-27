subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

!> Called to update q by solving source term equation 
!! $q_t = \psi(q)$ over time dt starting at time t.
!!
!! Integrate vertical velocity
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: i,j

    !write(6,*) '+++ in src2', maxval(q(6,:,:))

    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            ! integrate vertical velocity to get displacement:
            q(6,i,j) = q(6,i,j) + dt*q(5,i,j)
            enddo
        enddo

end subroutine src2
