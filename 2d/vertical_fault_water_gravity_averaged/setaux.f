c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # variable coefficient acoustics
c     #  aux(1,i,j) = density rho in (i,j) cell
c     #  aux(2,i,j) = lambda in (i,j) cell
c     #  aux(3,i,j) = mu in (i,j) cell
c     #  aux(4,i,j) = cp in (i,j) cell
c     #  aux(5,i,j) = cs in (i,j) cell
c     #  aux(6,i,j) = slip across fault
c     #  aux(7,i,j) = stretch in x for ABL
c     #  aux(8,i,j) = stretch in y for ABL
c
c     # Piecewise constant medium
c     # Material parameters are set in setprob.f

c
c     
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      common /comaux/ rho1,amu1,alam1,rho2,amu2,alam2

      real(kind=8) :: ABLdepth, ABLxpos(2), ABLypos
      common /ablparam/ ABLdepth, ABLxpos, ABLypos

      pi2 = 2.d0*datan(1.d0)

      aux(6,:,:) = 0.d0   ! set in b4step2 each time step
      aux(7,:,:) = 1.d0   ! reset below if ABL
      aux(8,:,:) = 1.d0   ! reset below if ABL

      do 30 j=1-mbc,my+mbc
       ycell = ylower + (j-0.5d0)*dy
       do 20 i=1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
          xl = xlower + (i-1.0d0)*dx
          yl = ylower + (j-1.0d0)*dy

          if (.true.) then
!$OMP         CRITICAL (cellave_fss)
              call cellave(xl,yl,dx,dy,w1)
!$OMP         END CRITICAL (cellave_fss)
           else
              w1 = 0.d0  ! no water, all rock
           endif

          w2 = 1.d0 - w1

          aux(1,i,j) = w1*rho1 + w2*rho2
          !aux(2,i,j) = w1*alam1 + w2*alam2
          aux(2,i,j) = 1.d0/(w1/alam1 + w2/alam2)
          aux(3,i,j) = w1*amu1 + w2*amu2

!         if (amu1*amu2 > 0.d0) then
!             aux(3,i,j) = 1.d0/(w1/amu1 + w2/amu2)
!           else if (w1==0.d0) then
!             aux(3,i,j) = amu2
!           else if (w2==0.d0) then
!             aux(3,i,j) = amu1
!           else 
!             aux(3,i,j) = 0.d0
!           endif

          !bulk1      = alam1 + 2.d0*amu1
          !bulk2      = alam2 + 2.d0*amu2
          !bulk       = 1.d0/(w1/bulk1 + w2/bulk2)
          bulk       = aux(2,i,j) + 2.d0*aux(3,i,j)
          aux(4,i,j) = dsqrt(bulk/aux(1,i,j))
          aux(5,i,j) = dsqrt(aux(3,i,j)/aux(1,i,j))

          ! set absorbing layer factor in x direction
          if (ABLdepth > 1.d-10) then
            if (xcell .le. ABLxpos(1)) then
              aux(7,i,j) = 1.d0/(1.d0 + dtan(pi2*(ABLxpos(1) 
     &                      - xcell)/ABLdepth)**2)
              aux(7,i,j) = exp(-2.d0*(ABLxpos(1)-xcell)/ABLdepth)
            elseif (xcell .ge. ABLxpos(2)) then
              aux(7,i,j) = 1.d0/(1.d0 + dtan(pi2*(xcell 
     &                      - ABLxpos(2))/ABLdepth)**2)
              aux(7,i,j) = exp(-2.d0*(xcell-ABLxpos(2))/ABLdepth)
            else
              aux(7,i,j) = 1.d0
            end if
            if (aux(7,i,j) > 1.d0) write(6,*) '+++ aux7 = ',aux(7,i,j)
          end if

          ! set absorbing layer factor in y direction
          if (ABLdepth > 1.d-10) then
            if (ycell .le. ABLypos) then
              aux(8,i,j) = 1.d0/(1.d0 + dtan(pi2*(ABLypos 
     &                     - ycell)/ABLdepth)**2)
              aux(8,i,j) = exp(-2.d0*(ABLypos-ycell)/ABLdepth)
            else
              aux(8,i,j) = 1.d0
            end if
            if (aux(8,i,j) > 1.d0) write(6,*) '+++ aux8 = ',aux(8,i,j)
          end if
   20     continue
   30    continue

       return
       end

