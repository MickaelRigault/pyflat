#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" This module defines the basics objects used for the star flat analysis"""

import warnings
import numpy     as np

# - Astropy
from astropy import units

# - Astrobject
from astrobject.utils.tools import kwargs_update
from astrobject.baseobject  import BaseObject
from astrobject.collections.photospatial import PhotoMapCollection


# - shapely
from shapely import geometry, vectorized


def get_flatfielder( photomaps, wcs_extension=None, **kwargs ):
    """ Get a PhotoMap collection

    Parameters
    ----------
    photomaps: [list of photomaps or filenames of]
        The photomaps of roughly a same area that you want to combine into a collection
        if a list of filename is given, get_sepobject() will be used to read them.
        This could be e.g. sep or sextractor outputs

    wcs_extension: [int/None] -optional-
        If the given data are fitfiles, you can specify the extension of the header
        containing the wcs solution if any. Leave this to None if no wcs solution have
        to be loaded.

    **kwargs goes to the PhotoMapCollection __init__ (catalogue etc.)
    
    Returns
    -------
    PhotoMapCollection
    """
    from astrobject import get_sepobject
    
    if not hasattr(photomaps, "__iter__"):
        raise TypeError("photomaps must be a list/array")

    # -- loading data
    if type(photomaps[0]) == str:
        photomaps = [get_sepobject(filename, wcs_extension=wcs_extension, **kwargs) for filename in photomaps]
        if wcs_extension is not None:
            [p._derive_radec_() for p in photomaps if p.has_wcs()]

    return FlatFielder(photomaps, **kwargs)

    
#######################################
#                                     #
#                                     #
#   Flat Fielder                      #
#                                     #
#                                     #
#######################################

class FlatFielder( PhotoMapCollection ):
    """ Collection gathering flatfielding observation """


    # ==================== #
    #  Main Methods        #
    # ==================== #
    # --------- #
    #  GETTER   #
    # --------- #
    def getcat_residual(self, param, catindex):
        """ get the residual parameter resulting from the
        difference between the parameter value in a given photomaps in
        comparison the its mean value accross the photomaps.
        
        Returns
        -------
        Arrays [NxM] where N is not number of photomaps and M the size of catindex
        """
        if hasattr(param, "__iter__"):
            raise TypeError("param must be a string not an array")
        if type(param) != str:
            raise TypeError("param must be a string (name of the parameter)")
        
        value = np.asarray(self.getcat(param, catindex))
        return value- np.nanmean(value, axis=0)
        
    # --------- #
    #  Project  #
    # --------- #
    def projec_to_grid(self):
        pass

    # --------- #
    #  PLOTTER  #
    # --------- #
    def show(self, ax=None, show=True, savefile=None,
             shownstep=1,
             cmap=None, catindexes=None, prop_catindex={}, **kwargs):
        """
        
        Parameters
        ----------

        shownstep: [int] -optional-
            Catalogue entries skipped from the plot ([::shownstep]).
            1 means no skipping
        
        **kwargs any mpl.scatter parameter
        
        Returns
        -------
        Void
        """
        import matplotlib.pyplot as mpl
        from astrobject.utils.mpladdon import figout

        #  Axis Definition
        # -----------------
        if ax is None:
            fig = mpl.figure(figsize=[5,8])
            ax  = fig.add_subplot(111)
        elif not hasattr(ax,"plot"):
            raise TypeError("The given ax is not a matplotlib axes")
        else:
            fig = ax.figure
            
        #  Definition
        # -----------------            
        prop= kwargs_update( dict(alpha = 0.7, mew=0, marker="o",ms=6, ls="None"),
                            **kwargs )
        if cmap is None:
            cmap = mpl.cm.viridis
            
        #  Data
        # -----------------
        if catindexes is None:
            if len(prop_catindex) ==0:
                warnings.warn("Default catatalog indexes loaded")
            catindexes = self.get_catindexes(**prop_catindex)
            
        x = self.getcat("x",catindexes[::int(shownstep)] if int(shownstep)>1 else catindexes )
        y = self.getcat("y",catindexes[::int(shownstep)] if int(shownstep)>1 else catindexes )
        
        [ax.plot(x[i],y[i], mfc=cmap(np.float(i)/(self.nsources-1)), **prop)
         for i in range(self.nsources)]
        
        #  Output
        # -----------------
        fig.figout(savefile=savefile, show=show)

    


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
    DERIVED_PROPERTIES = [] 

    def show(self, ax=None, cval=None, ec="0.5",
             savefile=None, show=True,
             vmin=None, vmax=None, cmap=None,
             cax=None,clabel="",cbarprop={}, scaleax=True,
             **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        from astrobject.utils.shape    import draw_polygon
        from astrobject.utils.mpladdon import insert_ax, colorbar, figout
        
        #  Axis Setting
        xmin,ymin, xmax, ymax = self.grid.bounds
        if ax is None:
            heigth = ymax-ymin
            width = xmax-xmin
            fig = mpl.figure(figsize=[5,5*float(heigth)/width])
            ax  = fig.add_subplot(111)
        else:
            fig = ax.figure
            
        # Colors
        if cval is not None:
            if len(cval) != self.npoly:
                raise ValueError("the size of cval (%d) do not corresponds to the number of polygons (%d)"%(len(cval),self.npoly))
            
            if vmin is None: vmin = np.percentile(cval[cval==cval], 5)
            if vmax is None: vmax = np.percentile(cval[cval==cval], 95)
            if cmap is None: cmap = mpl.cm.viridis
            cval_ = (cval-vmin)/(vmax-vmin)
            fcs = [cmap(c) if c==c else mpl.cm.binary(0.5,0.2) for c in cval_]
        else:
            fcs = ["None"]*self.npoly

        # Patches
        pol = [ draw_polygon(ax, p_, fc=fcs[i], ec=ec, **kwargs)
               for i,p_ in enumerate(self.grid)]
        
        # ColorBar
        if cax is None:
            cax = ax.insert_ax("right", shrunk=0.93,space=.0,axspace=0.03)
        if cax is not False:
            cax.colorbar(cmap,vmin=vmin,vmax=vmax,label=clabel,**cbarprop)

        # Scaling
        if scaleax:
            ax.set_xlim(xmin,xmax )
            ax.set_ylim(ymin,ymax )
            
        # Output
        fig.figout(savefile=savefile, show=show)

        
    def set_grid(self, multipolygon):
        """ """
        if geometry.MultiPolygon not in multipolygon.__class__.__mro__:
            raise TypeError("This input grid must be a shapely's MultiPolygon ")
        
        self._properties["grid"] = multipolygon
    
    def get_maskin(self, x, y):
        """ list of boolean mask saying if the points are wihtin
        the grid. Ordering following that of the grid (self.grid)
        """
        return [vectorized.contains(s,x,y) for s in self.grid]
    
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
