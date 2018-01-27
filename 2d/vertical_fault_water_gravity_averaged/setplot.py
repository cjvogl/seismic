
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
import sys
sys.path.insert(0,'..')

csig = 100
cdivcurl = 1e-6
xlimits = [-200e3, 200e3]

h0 = 4500.  # ocean depth
cw = 1500.  # speed of sound in water
g = 9.81
rho = 1025. # density of water
rhog = rho*g
period = 4 * h0 / cw # time for signal transit depth 4 times
ybot = -h0

x0_slope = 50e3
x0_shelf = 100e3
z0_ocean = -h0
#z0_shelf = z0_ocean #flat
z0_shelf = -960.0

def zbottom(x):
    #ybottom = -10000. + 6000*(x+40e3)/80e3

    slope_of_slope = (z0_ocean - z0_shelf) / (x0_slope - x0_shelf)
    ybottom = z0_ocean
    ybottom = np.where(x<x0_slope, z0_ocean, \
                    z0_ocean + slope_of_slope*(x-x0_slope))
    ybottom = np.where(x<x0_shelf, ybottom, z0_shelf)
    return ybottom

def plot_interfaces(current_data):
    from pylab import linspace, plot
    x = linspace(-200e3, 200e3, 101)
    plot(x, zbottom(x), 'k')

#--- Okada:

import setfault
reload(setfault)
fault = setfault.make_fault()
times = [1000.]  # after fault motion complete
x = np.linspace(-2,2,1000)
y = np.array([0,1])
dtopo = fault.create_dtopography(x,y,times)


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'
    

    def sigmatr(current_data):
        q = current_data.q
        return q[0,:,:] + q[1,:,:]

    def pressure(current_data):
        """ pressure in units of meters of water"""
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:])/(2.*rhog)

    def div(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            d = zeros(u.shape)
            return d
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        ux = (u[I+1,:][:,J] - u[I-1,:][:,J]) / (2*dx)
        vy = (v[:,J+1][I,:] - v[:,J-1][I,:]) / (2*dy)
        dint = ux + vy
        
        #zx = zeros((mx-2,1))
        #zy = zeros((1,my))
        #d = vstack((zy, hstack((zx, ux+vy, zx)), zy))
        
        d0 = dint[:,0]
        d1 = dint[:,-1]
        d2 = vstack((d0, dint.T, d1)).T
        d0 = d2[0,:]
        d1 = d2[-1,:]
        d = vstack((d0,d2,d1))      
        return d

    def curl(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            c = zeros(u.shape)
            return c
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        vx = (v[I+1,:][:,J] - v[I-1,:][:,J]) / (2*dx)
        uy = (u[:,J+1][I,:] - u[:,J-1][I,:]) / (2*dy)
        cint = vx - uy

        c0 = cint[:,0]
        c1 = cint[:,-1]
        c2 = vstack((c0, cint.T, c1)).T
        c0 = c2[0,:]
        c1 = c2[-1,:]
        c = vstack((c0,c2,c1))      
        return c

    # Figure for trace(sigma) and sigma_12 side by side
    plotfigure = plotdata.new_plotfigure(name='P and S waves', figno=11)
    plotfigure.kwargs = {'figsize':(13,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.6,.35,.35])' # 'subplot(221)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'pressure'
    #plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = pressure
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -csig
    plotitem.pcolor_cmax = csig
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.5,.6,.45,.35])' # 'subplot(222)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'delta' #'sigma_12'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 5 #2
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.5 #-csig
    plotitem.pcolor_cmax =  0.5 #csig
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.1,.35,.35])' # 'subplot(223)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'div(u)'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = div
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -cdivcurl
    plotitem.pcolor_cmax = cdivcurl
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]


    # Figure for u-velocity:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.show = False
    plotaxes.axescmd = 'axes([.5,.1,.45,.35])' # 'subplot(224)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'u-velocity'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 3
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.002
    plotitem.pcolor_cmax = 0.002
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]


    # Figure for curl:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.show = False
    plotaxes.axescmd = 'axes([.5,.1,.45,.35])' # 'subplot(224)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'curl(u)'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = curl
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -cdivcurl
    plotitem.pcolor_cmax = cdivcurl
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]



    # Figure for surface displacement
    plotfigure = plotdata.new_plotfigure(name='surface displacement', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [0,2]
    plotaxes.ylimits = [-0.8,0.8]
    plotaxes.title = 'surface/bottom displacement'
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.show = False

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        global surf_min, surf_max
        from pylab import find,ravel,nan
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= 0.) & (y > -dy))
        x_slice = ravel(x)[ij]

        # extrapolate displacement from cell centers to top surface:
        eta_slice_1 = ravel(q[5,:,:])[ij]
        ij = find((y <= -dy) & (y > -2*dy))
        eta_slice_2 = ravel(q[5,:,:])[ij]
        eta_slice = 1.5*eta_slice_1 - 0.5*eta_slice_2

        #if current_data.level < 3:
            #eta_slice = nan*eta_slice
        #if len(eta_slice) > 0:
            #print('min,max = ',eta_slice.min(),eta_slice.max())
        try:
            surf_min = min(surf_min,eta_slice.min())
            surf_max = max(surf_max,eta_slice.max())
        except:
            pass
        if current_data.level != 2:
            eta_slice *= nan
        return x_slice, eta_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = 'b-'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_data_show = [1]  

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.show = False

    def xsec_bottom_delta(current_data):
        from pylab import find,ravel,nan
        
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= ybot+dy) & (y > ybot))
        x_slice = ravel(x)[ij]

        # bottom displacement:
        delta_slice = ravel(q[5,:,:])[ij]
        if current_data.level != 2:
            delta_slice *= nan
        return x_slice, delta_slice

    plotitem.map_2d_to_1d = xsec_bottom_delta
    plotitem.plotstyle = 'g-'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_data_show = [1]  

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.show = False

    def xsec_bottom_pressure(current_data):
        from pylab import find,ravel,nan
        
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= ybot+dy) & (y > ybot))
        x_slice = ravel(x)[ij]

        # bottom pressure in units of meters of water, corrected for delta:
        p_slice = -(ravel(q[0,:,:])[ij] + ravel(q[1,:,:])[ij]) / (2*rhog) \
                    - ravel(q[5,:,:])[ij]
        if current_data.level != 2:
            p_slice *= nan
        return x_slice, p_slice





    plotitem.map_2d_to_1d = xsec_bottom_pressure
    plotitem.plotstyle = 'r-'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_data_show = [1]  

    def plot_bottom_pressure(current_data):
        from numpy import linspace, ones
        from pylab import plot, subplot, xlim, ylim, title, legend
        from clawpack.visclaw import gridtools

        framesoln = current_data.framesoln
        xout = linspace(xlimits[0], xlimits[1], 1000)
        #eta_swe, u_swe = swe_solution(xout, framesoln.t)
        yout = zbottom(xout) + 10.
        delta_bottom = gridtools.grid_output(framesoln, 5, xout, yout)
        pressure = -gridtools.grid_output(framesoln, 0, xout, yout)
        pdiff = pressure  / rhog - delta_bottom
        #plot(xout, pdiff, 'r')
        plot(111e3*dtopo.X[0,:],dtopo.dZ[-1,0,:],'m')
        plot(xout, delta_bottom, 'g')
        #import pdb; pdb.set_trace()


    def fixup(current_data):
        from pylab import grid
        from clawpack.visclaw import legend_tools

        plot_bottom_pressure(current_data)
        grid(True)
        #labels = ['surface displacement','bottom displacement',
        #          'bottom pressure diff']
        #colors = ['b','g','r']
        labels = ['surface displacement','bottom displacement',
                  'Okada final']
        colors = ['b','g','m']
        linestyles = '-'
        markers = ''
        legend_tools.add_legend(labels,colors,linestyles,markers,
                loc='upper left')
        
    plotaxes.afteraxes = fixup

    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'Grid patches'
    plotaxes.scaled = False

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [1,0,0]
    plotitem.amr_patchedges_show = [1]


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    def aa_gauge(current_data):
        from pylab import grid, tight_layout
        grid(True)
        tight_layout()

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Horizontal velocity'
    plotaxes.afteraxes = aa_gauge

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Vertical velocity'
    plotaxes.afteraxes = aa_gauge

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 4
    plotitem.plotstyle = 'b-'


    plotfigure = plotdata.new_plotfigure(name='gauge plot sigma', figno=301, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    def pressure_diff(current_data):
        q = current_data.q
        p = -(q[0,:]+q[1,:]) / (2*rhog)
        delta = q[5,:]
        return p - delta

    def filtered_pressure_diff(current_data):
        from pylab import find,len,sum,zeros,hstack
        q = current_data.q
        import pdb; pdb.set_trace()
        p = -(q[0,:]+q[1,:]) / (2*rhog)
        delta = q[5,:]
        d = p - delta
        t = q.t
        fpd = []
        for k1 in range(len(t)):
            t1 = t[k1]
            try:
                k2 = find(t > t1+period).min()
            except:
                break
            if k1==0:
                kmid = int((k1+k2)/2.)
            fpd.append(sum(d[k1:k2])/(k2-k1+1))
        fpd = hstack((zeros(kmid),fpd))
        kz = len(t) - len(fpd)
        if kz>0:
            fpd = hstack((fpd,zeros(kz)))
        return fpd

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'pressure difference'
    plotaxes.afteraxes = aa_gauge

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = pressure_diff
    plotitem.plotstyle = 'b-'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = filtered_pressure_diff
    plotitem.plotstyle = 'k-'



    def displacement(current_data):
        q = current_data.q
        delta = q[5,:]
        return delta

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'displacement'
    plotaxes.afteraxes = aa_gauge

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = displacement
    plotitem.plotstyle = 'b-'



    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True

    return plotdata

    
