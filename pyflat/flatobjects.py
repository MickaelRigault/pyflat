#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" This module defines the basics objects used for the star flat analysis"""

import warnings
import numpy     as np
from scipy import stats
import matplotlib.pyplot as mpl

# - Astropy
from astropy import units

# - Astrobject
from astrobject.utils.tools import kwargs_update
from astrobject.baseobject  import BaseObject
from astrobject.collections.photospatial import PhotoMapCollection


# - shapely
from shapely import geometry, vectorized

# - deco
#from deco import concurrent, synchronized

# - Local


    
#######################################
#                                     #
#                                     #
#   Flat Fielder                      #
#                                     #
#                                     #
#######################################

class FlatFielder( PhotoMapCollection ):
    """ Collection gathering flatfielding observation """

    DERIVED_PROPERTIES = ["trashed", "timecorr"]
    
    # ==================== #
    #  Main Methods        #
    # ==================== #
    # --------- #
    #  IO       #
    # --------- #
    # --------- #
    #  MANIP    #
    # --------- #
    def clean(self, catindex,verbose=True,
              min_stars=100, max_magres=0.005,
              savefile=None, show=False):
        """ """
        if min_stars is None and max_magres is None:
            warnings.warn("Nothing to clean")
            return
        
        
        oldid = self.list_id.copy()

        if min_stars is not None:
            delta_mag = np.asarray(self.getcat_residual("mag_auto", catindex))
            rm_ = [self.remove(id_) for id_ in self.list_id[np.sum(1-np.isnan(delta_mag.T), axis=0)<100]]
            if verbose:
                print "%d removed for they have less than 100 stars: "%len(rm_) +", ".join([id_ for id_ in oldid if id_ not in self.list_id])

                
        if max_magres is not None:
            delta_mag_corr, corr = self.get_global_magtrend("mjd", catindex, max_magres=max_magres,
                                                            savefile=savefile, show=show)
            rm_ = [self.remove(id_) for id_ in self.list_id[np.abs(np.nanmedian(delta_mag_corr.T,axis=0))>max_magres]]
            if verbose:
                print "%d removed for they large residual"%len(rm_)



    def load_nightly_evolution(self, catindex=None, magkey="mag_auto",
                               show=True):
        """ Return the night evolution based on stellar magnitudes.
        
        Parameters
        ----------
        catindex: [N*M mask]
        """
        print "Not Done YET"
        
        #if catindex is None:
       #     catindex = self.get_catindexes(1, isolated_only=True, stars_only=True)
        # - get the data you need
       # magres  = np.asarray(self.getcat_residual(magkey, catindex, nanflagged=True))
       # mjdres  = np.asarray(self.getcat_residual("mjd", catindex, clipping=False))
       # if show:
       #     fig = mpl.add_subplot()
        
        
    
        
        
        




    



                
    def get_global_magtrend(self, key, catindex, magkey ="mag_auto",
                            max_magres=None, savefile=None, show=False):
        """ """
        delta_mag = np.asarray(self.getcat_residual(magkey, catindex))
        corrv       = np.asarray(self.getcat_residual(key, catindex))
        from .analysis import get_correction
        if savefile is not None or show:
            fig = mpl.figure(figsize=[8,6])
            ax  = fig.add_subplot(111)
            ax.set_xlabel(r"$\mathrm{\Delta\ %s}$"%key, fontsize="xx-large")
            ax.set_ylabel(r"$\mathrm{global\ magnitude\ residual}$", fontsize="xx-large")
        else:
            ax = None
            
        return get_correction(np.nanmean(corrv.T, axis=0), np.nanmedian(delta_mag.T,axis=0),
                                  ax=ax, show=show, highlight_extrem=max_magres), delta_mag
        

    
    def set_magcorr(self, key, metakey, catindex, magkey="mag_auto", savefile=None, show=False,
                    maxres=None):
        """ correct the magnitude from its global trend as a function of `key`
        This will be saved as `metakey`
        """
        delta_mag_corr, corr = self.get_global_magtrend(key, catindex,magkey=magkey,
                                                        savefile=savefile, show=show,  
                                                        max_magres=maxres)
        self.set_meta(metakey, delta_mag_corr, catindex)
        
        
    # --------- #
    #  GETTER   #
    # --------- #
    def getcat_residual(self, param, catindex, nanflagged=False, usemedian=True,
                        clipping = False,  sigma_low=4, sigma_high=4):
        """ get the residual parameter resulting from the
        difference between the parameter value in a given photomaps in
        comparison the its mean value accross the photomaps.
        
        Returns
        -------
        Arrays [NxM] where N is not number of photomaps and M the size of catindex
        """
        from .utils import get_residual
        if hasattr(param, "__iter__"):
            raise TypeError("param must be a string not an array")
        if type(param) != str:
            raise TypeError("param must be a string (name of the parameter)")
        
        # - value
        value = np.asarray( self.getcat(param, catindex) )
        if nanflagged:
            value[np.asarray(self.getcat("flags", catindex))] = np.NaN
        if not clipping:
            sigma_low = sigma_high = 1e5
        return get_residual(value, sigma_low, sigma_high, usemedian=usemedian)


    # --------- #
    #  SETTED   #
    # --------- #
    def set_meta(self, key, value, catindex=None):
        """ Set the Meta to each PhotoMaps for the given catindex"""
        if not hasattr(value,"__iter__"):
            value = [value] *len(self.list_id)
        elif len(value) != len(self.list_id):
            raise TypeError("the given value must have a same size as list_id")
        
        if catindex is None:
            [self.photomaps[id_].set_meta(key,v_) for id_,v in zip(self.list_id,value)]
        else:
            for pmap,v_ in zip(self.photomaps.values(), value):
                idx = pmap.catindex_to_index(catindex, cleanindex=True)
                valuein = np.ones(len(pmap.list_id))*np.NaN
                valuein[idx] = v_
                pmap.set_meta(key, valuein)
        
    
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

    
    @property
    def list_id(self):
        """ sorted as int """
        return np.asarray(np.sort(np.asarray( super(FlatFielder, self).list_id, dtype="int")), dtype=str)




        
#######################################
#                                     #
#                                     #
#   Grids                             #
#                                     #
#                                     #
#######################################
