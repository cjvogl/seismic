c
c
c
c     =================================================
      function fdisc(x,y)
c     =================================================

      implicit double precision (a-h,o-z)

c
c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the 
c     # left of the curve and positive to the right

c     # fdisc is positive in the lower layer, negative in upper

      !fdisc = -10000. + 6000.*(x+200.d3)/400.d3  - y
      !fdisc = -4500.d0 - y

      ! flat ocean with slope up to flat shelf, piecewise linear:
      x0_slope = 50d3
      x0_shelf = 100d3
      z0_ocean = -4500.d0
      !z0_shelf = z0_ocean !flat
      z0_shelf = -960.d0
      ybottom = -4500.d0
      slope_of_slope = (z0_ocean - z0_shelf) / (x0_slope - x0_shelf)
      ybottom = z0_ocean
      if (x > x0_slope) 
     &   ybottom = z0_ocean + slope_of_slope*(x-x0_slope)
      if (x > x0_shelf) ybottom = z0_shelf

      !ytop = 0.d0
      !yp = ybottom + (yc-yc_lower)/(yc_upper-yc_lower) * (ytop - ybottom)
      
      fdisc = ybottom - y
c
      return
      end

