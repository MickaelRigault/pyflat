#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" This module defines the basics objects used for the star flat analysis"""

import warnings
import numpy     as np

from shapely import geometry, vectorized
# - Astrobject
from astrobject.baseobject  import BaseObject
from astrobject.utils       import shape

def get_voronoy_grid(x, y, edges, npergrid=None):
    """ """
    flagnan = np.isnan(x) * np.isnan(y)
    index_used = np.arange(len(x))[~flagnan]

    # - Regular Voronoy
    if npergrid is None or npergrid<=1:
        mpoly_voronoi = shape.get_voronoy_multipolygon(x[index_used],y[index_used], edges)
    else:
        from scipy.stats import gaussian_kde
        xy = np.vstack([x[index_used],y[index_used]])
        z = gaussian_kde(xy)
        xvor,yvor = z.resample(int(len(x[index_used])/float(npergrid)))
        mpoly_voronoi = shape.get_voronoy_multipolygon(xvor,yvor, edges)
    
    return PatchGrid(mpoly_voronoi)


def get_square_grid(rect, ngrid):
    """
    Parameters:
    -----------
    rect: [4 floats/ints]
        The total grid area: [xmin, ymin, xmax, ymax]
        
    ngrid: [2 ints]
        x and y binning.

    Returns
    -------
    PatchGrid
    """
    baserect= np.asarray([[0,0],[0,1],[1,1], [1,0]])
    # - prop
    xmin,ymin, xmax, ymax = rect
    width  = float(xmax-xmin)
    heigth = float(ymax-ymin)
    xbin, ybin = ngrid 
    # - single rect
    scaling = np.asarray([width/xbin,heigth/ybin])
    sourcing = np.asarray([xmin, ymin])

    polys = []
    rect_p_source = baserect * scaling
    for i in range(xbin):
        xoffset = i*scaling[0]
        for j in range(ybin):
            yoffset = j*scaling[1]
            rect_p = np.asarray(rect_p_source + [[xoffset,yoffset]]*4)
            polys.append( geometry.Polygon(rect_p))
        
    mpoly = geometry.MultiPolygon(polys)
    
    return PatchGrid(mpoly)

#######################################
#                                     #
#                                     #
#   Grids                             #
#                                     #
#                                     #
#######################################
class PatchGrid( BaseObject ):
    """ """
    PROPERTIES         = ["grid"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["masking"] 

    def __init__(self, mpoly=None, empty=False):
        """ Initialize the Grid """
        self.__build__()
        if empty:
            return
        if mpoly is not None:
            self.set_grid(mpoly)

    # ----------- #
    #  MASKING    #
    # ----------- #
    def set_data_masking(self, x, y):
        """ attach to this instance the masking information of the given data """
        self._derived_properties["masking"] = self.get_cell_masking(x, y)
        
    def get_cell_masking(self, x, y):
        """ Returns a boolean masking for each polygon on weither or not
        it contains the given x, y points """
        return np.asarray([vectorized.contains(s_, x, y) for s_ in self.grid])

    def project(self, value):
        """ Returns the dictionary of value associated with the ith polygon """
        if not self.has_flagin():
            raise AttributeError("The masking has not been set. run self.set_data_masking() ")
        value = np.asarray(value)
        return [value[self.flagin[i]] for i in range(self.npoly)]
    

    def show(self, cval=None, ax=None, ec="0.5",
             savefile=None, show=True,
             cax=None, clabel="",cfontsize="x-large",
             cbarprop={}, scaleax=True,
             **kwargs):
        """ """
        from .utils import polygonplot
        import matplotlib.pyplot as mpl
        #  Axis Setting
        xmin,ymin, xmax, ymax = self.grid.bounds
        if ax is None:
            heigth = ymax-ymin
            width = xmax-xmin
            fig = mpl.figure(figsize=[5,5*float(heigth)/width])
            ax  = fig.add_subplot(111)
            
        return ax.polygonplot(self.grid, cval=cval, savefile=savefile, show=show,
                              ec=ec, cax=cax, clabel=clabel, cfontsize=cfontsize,
                              scaleax=scaleax, cbarprop=cbarprop,  **kwargs)



    def show_mag_pull(self,data, fig= None, **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        from scipy import stats
        
        if fig is None:
            fig = mpl.figure(figsize=[10,10])
    
        ax   = fig.add_axes([0.1,0.1,0.37,0.8])
        axh  = fig.add_axes([0.56,0.5,0.37,0.25])

        self.show(cval = data, ax=ax, show=False, savefile=None,
                  clabel=r"$\mathrm{data/error}$",
                  **kwargs)

        # -- Histo
        prop = dict(histtype="step", fill=True, fc=mpl.cm.binary(0.6,0.4),
                    ec=mpl.cm.binary(0.5,0.9), 
                    lw=0, bins=20, normed=True, range=[-8,8])
        # Data
        axh.hist(data, **prop)
        mean_, std_ = np.nanmean(data),np.nanstd(data)
        # Expected
        axht = axh#.twinx()
        xx = np.linspace(-7,7,1000)
        
        axht.plot( xx, stats.norm.pdf(xx,loc=0, scale=1), color="k", lw=2, label=r"$\mathrm{Expected}$")
        axht.plot( xx, stats.norm.pdf(xx,loc=mean_, scale=std_), color="k", ls="--", label=r"$\mathrm{Observed}$")
        
        axht.set_title(r"$\mathrm{Pull\ Distribution:\ x_0=%.2f\ \sigma=%.2f}$"%(mean_, std_), 
                    fontsize="x-large")
        axh.set_xlabel(r"$\mathrm{data/error}$",fontsize="xx-large")
        axh.set_ylim(0, axh.get_ylim()[1]*1.2)
        axh.legend(loc="upper right", ncol=1, frameon=False)
        
        _ = [ax_.set_yticks([]) for ax_ in [axht,axh]]
        
    
    def set_grid(self, multipolygon):
        """ """
        if geometry.MultiPolygon not in multipolygon.__class__.__mro__:
            raise TypeError("This input grid must be a shapely's MultiPolygon ")
        
        self._properties["grid"] = multipolygon
    
    
    # ================== #
    #  Properties        #
    # ================== #
    @property
    def grid(self):
        """ """
        return self._properties["grid"]
    
    def has_grid(self):
        """ Does this method has a grid loaded? """
        return self.grid is not None
    
    @property
    def npoly(self):
        """ number of polygon in the grid """
        if not self.has_grid():
            return None
        return len(self.grid)
    @property
    def flagin(self):
        """ boolean flag indicating which point is in which grid """
        return self._derived_properties["masking"]

    def has_flagin(self):
        """ Test if the masking has been set """
        return self.flagin is not None
