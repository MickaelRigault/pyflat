#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Tools to analysis the data"""

import warnings
import numpy     as np

import matplotlib.pyplot as mpl
from scipy.special import orthogonal
from scipy.optimize import curve_fit
from astrobject.utils.tools import kwargs_update
###############################
#                             #
# Fitting and Removing Trend  #
#                             #
###############################
def get_correction(x, y, degree=2, verbose=False,
                   savefile=None,show=False,ax = None,
                   highlight_extrem=None, 
                   **kwargs):
    """ This is a model for the evolution. Do y - what_is_returned
    
    """
    def first_order(x_, a, b):
        return a + b*x_ 
    def second_order(x_, a, b, c):
        return a + b*x_ + c*x_**2
    def third_order(x_, a, b, c, d):
        return a + b*x_ + c*x_**2 + d*x**3
    
    if degree >3:
        raise ValueError("only 1, 2 or 3 degrees done")
    if degree==1:
        polycurve = first_order
    elif degree==2:
        polycurve = second_order
    elif degree==3:
        polycurve = third_order
        
    # - clean the data
    flagnan = np.isnan(x) * np.isnan(y)
    xused,yused = x[~flagnan], y[~flagnan]
    popt, pcov = curve_fit(polycurve, xused, yused)
    if verbose:
        print popt
        
    if show or savefile is not None or ax is not None:
        from astrobject.utils.mpladdon import figout

        if ax is None:
            fig = mpl.figure(figsize=[8,5])
            ax  = fig.add_subplot(111)
        else:
            fig = ax.figure
            
        prop = kwargs_update( dict(ms=15, ls="None",mfc=mpl.cm.Blues(0.6,0.5),
                                   mec=mpl.cm.Blues(0.6), mew=2,zorder=4,marker="o"),
                                    **kwargs)
        
        if highlight_extrem is None:
            _ = ax.plot(xused,yused,  **prop)
        else:
            ymodel = polycurve(xused,*popt)
            flagin = np.abs(yused-ymodel)<=highlight_extrem
            
            
            new_prop = kwargs_update(prop,**dict(mfc=mpl.cm.bone(0.6,0.5),mew=0,
                                     mec=mpl.cm.bone(0.6),zorder=prop["zorder"]-1))
            
            _ = ax.plot(xused[flagin],yused[flagin], 
                           **prop)
            
            _ = ax.plot(xused[~flagin],yused[~flagin], 
                           **new_prop)
            

            
        xx = np.linspace(xused.min(),xused.max(), 100)    
        ax.plot(xx,polycurve(xx,*popt) , marker="None", ls="-",
                scalex=False, scaley=False,color="k", lw=2, zorder=prop["zorder"]+1)
        
        fig.figout(savefile=savefile, show=show)
        
    return polycurve(x,*popt)

    



