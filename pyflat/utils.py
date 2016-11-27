#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" ToolBox for the module """

import numpy as np
import warnings
from scipy import stats

from astrobject.utils.decorators import make_method
import matplotlib.pyplot as mpl


# ======================== #
#   Array Manipulation     #
# ======================== #

def get_residual(value, sigma_low=4, sigma_high=4,
                 usemedian=True):
    """ Returns the residual: from flatobject's star/time structure
                 the average/median magnitude of a given star is removed
                 per star
    Returns
    -------
    Array (same structure as value)
    """
    if usemedian:
        return value - np.asarray([ np.nanmedian( stats.sigmaclip( v[~np.isnan(v)],sigma_low,sigma_high)[0]) for v in value.T])
    return value - np.asarray([ np.nanmean( stats.sigmaclip( v[~np.isnan(v)],sigma_low,sigma_high)[0]) for v in value.T])

# ======================== #
#   MPL addon              #
# ======================== #
@make_method(mpl.Axes)
def polygonplot(ax, polygones, cval=None, cmap=mpl.cm.viridis,
                ec="0.5",savefile=None, show=True,
                vmin=None, vmax=None, cax=None, clabel="", cfontsize="x-large",
                cbarprop={}, scaleax=True,**kwargs):
    
    """ Displays shapely's polygons / multipolygon on the axis """
    
    from shapely import geometry
    from astrobject.utils.shape    import draw_polygon
    from astrobject.utils.mpladdon import insert_ax, colorbar, figout
        
    # Input
    if geometry.multipolygon.MultiPolygon not in polygones.__class__.__mro__:
        # - This should be a multipolygon, mauybe you gave a list of polygon:
        if shapely.geometry.polygon.Polygon in polygones.__class__.__mro__:
            # Yes you did
            try:
                polygons = geometry.multipolygon.MultiPolygon(polygons)
            except:
                raise TypeError("Cannot convert polygons into MultiPolygon")
        else:
            TypeError("The given polygon should be a shapely's MultiPolygon or list of polygon")
    
    # Color
    if cval is not None:
        if len(cval) != len(polygones):
            raise ValueError("the size of cval (%d) do not corresponds to the number of polygons (%d)"%(len(cval), len(polygones)))
        cval = np.asarray(cval)
        if vmin is None: vmin = np.percentile(cval[cval==cval], 5)
        if vmax is None: vmax = np.percentile(cval[cval==cval], 95)
        if cmap is None: cmap = mpl.cm.viridis
        cval_ = (cval-vmin)/(vmax-vmin)
        fcs = [cmap(c) if c==c else "None" for c in cval_]
    else:
        fcs = ["None"] * len(polygones)

    # Patches
    pol = [ draw_polygon(ax, p_, fc=fcs[i], ec=ec, **kwargs)
            for i,p_ in enumerate(polygones)]
    # ColorBar
    if cax is None and cval is not None:
        cax = ax.insert_ax("right", shrunk=0.93,space=.0,axspace=0.01)
    else:
        cax = False
        
    if cax is not False:
        print cax
        cax.colorbar(cmap,vmin=vmin,vmax=vmax,
                     label=clabel, fontsize=cfontsize,
                     **cbarprop)

    # Scaling
    if scaleax:
        xmin,ymin, xmax, ymax = polygones.bounds
        ax.set_xlim(xmin,xmax )
        ax.set_ylim(ymin,ymax )
            
    # Output
    ax.figure.figout(savefile=savefile, show=show)
    
