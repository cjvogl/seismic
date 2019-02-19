
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from clawpack.clawutil.data import ClawData
from mapping import Mapping
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
from clawpack.geoclaw.data import LAT2METER
reload(dtopotools)

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    fault = dtopotools.Fault()
    fault.read(plotdata.outdir + '/fault.data')
    probdata = ClawData()
    probdata.read(plotdata.outdir + '/setprob.data',force=True)
    clawdata = ClawData()
    clawdata.read(plotdata.outdir + '/claw.data', force=True)

    mapping = Mapping(fault, probdata, clawdata)
    fault_width = mapping.fault_width
    xcenter = mapping.xcenter
    zcenter = mapping.zcenter
    xp1 = mapping.xp1
    xp2 = mapping.xp2
    zp1 = mapping.zp1
    zp2 = mapping.zp2

    xlimits = [xcenter-0.5*probdata.domain_width+probdata.abl_depth,
               xcenter+0.5*probdata.domain_width-probdata.abl_depth]
    zlimits = [-probdata.domain_depth+probdata.abl_depth, 0.0]
    xlimitsW = xlimits
    zlimitsW = [probdata.zlower_ocean, 0.0]
    abl_depth = probdata.abl_depth

    gaugedata = ClawData()
    gaugedata.read(plotdata.outdir + '/gauges.data',force=True)
    ngauges = gaugedata.ngauges/2
    xc = np.zeros(ngauges)
    for j in range(ngauges):
        g = plotdata.getgauge(j)
        xc[j] = g.location[0]
    fault.create_dtopography(xc/LAT2METER,np.array([0.]),[1.0],horiz_disp=True)

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'
    plotdata.parallel = True

    def plot_interfaces(current_data):
        from pylab import linspace, plot
        xl = (xp1,xp2)
        yl = (zp1,zp2)
        plot(xl,yl,'g',linewidth=2)
        xl = (xlimits[0],xlimits[1])
        yl = (0.0,0.0)
        plot(xl,yl,'b',linewidth=2)
        xl = (xlimits[0],mapping.xlower_slope,mapping.xlower_shelf,xlimits[1])
        yl = (mapping.zlower_ocean,mapping.zlower_ocean,mapping.zlower_shelf,mapping.zlower_shelf)
        plot(xl,yl,'k',linewidth=2)

    def after_ocean_pressure(current_data):
        from pylab import gca
        plot_interfaces(current_data)
        gca().set_ylabel('ocean pressure',fontsize=18)


    def plot_seasurface(current_data):
        from pylab import plot,zeros,gca, legend
        t = current_data.t

        zs = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(ngauges+gaugeno)
            for k in range(1,len(g.t)):
                if (g.t[k] > t and g.t[k-1] <= t):
                    dt = g.t[k] - g.t[k-1]
                    zs[gaugeno] = g.q[5,k]*(t-g.t[k-1])/dt +\
                                  g.q[5,k-1]*(g.t[k]-t)/dt

        ax = gca()
        plot(xc[:ngauges],zs,'-b',linewidth=2)
        ax.set_ylabel('ocean surface',fontsize=18)

    def plot_seafloor(current_data):
        from pylab import plot,zeros,gca, legend, colorbar
        t = current_data.t

        zf = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(gaugeno)
            for k in range(1,len(g.t)):
                if (g.t[k] > t and g.t[k-1] <= t):
                    dt = g.t[k] - g.t[k-1]
                    zf[gaugeno] = g.q[5,k]*(t-g.t[k-1])/dt +\
                                  g.q[5,k-1]*(g.t[k]-t)/dt

        ax = gca()
        plot(xc[:ngauges],zf,'-k')
        ax.set_ylabel('ocean floor',fontsize=18)

    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:])

    def downdip_vel(current_data):
        # return vel dot tau, where tau is tangent to fault
        tau_x = (xp2 - xp1)/fault_width
        tau_y = (zp2 - zp1)/fault_width
        u = current_data.q[3,:,:]
        v = current_data.q[4,:,:]
        return u*tau_x + v*tau_y

    def after_downdip_velocity(current_data):
        from pylab import gca, gcf
        plot_interfaces(current_data)
        ax = gca()
        ax.set_xlabel('m',fontsize=18)
        ax.set_ylabel('slip-direction velocity',fontsize=18)
        f = gcf()

    # Figure
    plotfigure = plotdata.new_plotfigure(name='results', figno=1)
    plotfigure.kwargs = {'figsize':(8,12), 'tight_layout': True}

    # Set axes for seasurface:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(411)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = [-0.3,0.5]
    plotaxes.title = ''
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_seasurface

    # Set axes for pressure in ocean:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(412)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = zlimitsW
    plotaxes.title = ''
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = after_ocean_pressure
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -1e5
    plotitem.pcolor_cmax = 1e5
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    # Set axes for seafloor:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(413)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = [-0.3,0.5]
    plotaxes.title = ''
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_seafloor

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(414)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = zlimits
    plotaxes.title = ''
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = after_downdip_velocity

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = downdip_vel
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.1
    plotitem.pcolor_cmax = 0.1
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

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
#    plotdata.parallel = True

    return plotdata
