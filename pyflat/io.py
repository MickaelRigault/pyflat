#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" This module contains the IO of the system """

import os
import numpy as np
import warnings

from astropy.io import fits
from astropy import table

import astrobject

DATAROOT = os.getenv("DATAPATH",None)
if DATAROOT is None:
    warnings.warn("YOU DATAPATH has not been defined. Some IO functionality won't work")



__all__ = ["get_filenames","get_sepmap","get_flatfielder"]


    
def get_filenames(ccd="*", night="*"):
    """ get the path of the files containing the given entries
    Parameters
    ----------
    ccd: [string]
        number of the ccd. *_c`ccd`* will be searched for, caution you might
        need to say e.g., ccd='01'

    night: [string]
        The filename should contain this name prior the ccd identification

    Returns
    -------
    list 
    """
    if DATAROOT is None:
        raise IOError("No environment variation 'DATAPATH' defined")
    
    from glob import glob
    return glob(DATAROOT+"PTF/ptf_starflat/catalogues/*%s*_c%s*"%(night,ccd))


def get_sepmap(filename, fluxkey="flux_auto", fluxkeyerr="fluxerr_auto",
               get_table=False, wcs_extension=None):
    """ """
    # Internal Tools
    # ---------------
    f = fits.open(filename)
    data,h2,h3 =  f[1].data,f[2].header,f[3].header
    
    def sex_to_sep_rename(t_):
        """ This define the baseline for PhotoMaps. Fast (< 1ms for ~4000 entries)"""
        keysswitch = {"ra":"x_world","dec":"y_world",
                      "x":"x_image","y":"y_image",
                      "a":"a_image","b":"b_image",
                      "a.err":"erra_image","b.err":"errb_image",
                      "theta":"theta_image","theta.err":"thetaerr_image",
                      "mjd":"obsmjd"}
        [t_.rename_column(v,k) for k,v in keysswitch.items() if v in t_.columns]

    hnames = ["CCDID","PTFFIELD","FILTER","FILTERID",
               "OBSMJD","SEEING","FWHMSEX","FBIAS","MSMAPCZP","AZIMUTH","ALTITUDE",
               "AIRMASS","OUTTEMP","OUTRELHU","GAIN","READNOI","DARKCUR",
               'AEXPTIME', 'APSFILT', 'APSCOL', 'APRMS', 'APBSRMS', 'APNZPRMS',
               'APNSTDI1', 'APNSTDIF', 'APCHI2', 'APDOF', 'APMEDJD', 'APPAR01', 'APPARE01',
               'APPAR02', 'APPARE02', 'APPAR03', 'APPARE03', 'APPAR04', 'APPARE04', 'APPAR05',
               'APPARE05', 'APPAR06', 'APPARE06', 'APPAR07', 'APPARE07', 'APPAR08', 'APPARE08',
               'APPAR09', 'APPARE09', 'APPAR10', 'APPARE10', 'APPAR11', 'APPARE11']
    
    def read_header(t_, h_, names):
        """ new entries based on header """
        nt_ = len(t_)
        t_.add_columns([table.Column([h_[name]]*nt_, name.lower()) 
                        for name in names if name in h_])
        
    # The Jobs
    # ---------------
    # Open the data
    # Convert that in table
    t = table.Table(data, names=[l.lower() for l in data.columns.names])
    # Add the missing keys
    read_header(t, h2, hnames)
    read_header(t, h3, hnames)
    # Rename keys
    sex_to_sep_rename(t)
    # Which is the flux?
    t.add_column(table.Column(t[fluxkey],"flux"))
    t.add_column(table.Column(t[fluxkeyerr]**2,"var"))

    if get_table:
        f.close()
        return t
    
    # - wcs if Needed
    wcs  = None if wcs_extension is None else astrobject.astrometry.wcs(header=f[wcs_extension].header)
    # The actual object
    s_ = astrobject.photospatial.SepObject()
    s_.create_from_table(t)
    if wcs is not None:
        s_.set_wcs(wcs)
        
    f.close()
    return s_
    
def get_flatfielder( filephotomaps, wcs_extension=None, **kwargs ):
    """ Get a PhotoMap collection

    Parameters
    ----------
    photomaps: [list of photomaps files]
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
    from .flatobjects import FlatFielder
    print "Flat fielder loading"
    if not hasattr(filephotomaps, "__iter__"):
        raise TypeError("photomaps must be a list/array. For Single file loading use get_sepmap")

    # -- loading data
    if type(filephotomaps[0]) == str:
        photomaps = [get_sepmap(filename, wcs_extension=wcs_extension, **kwargs) for filename in filephotomaps]
    else:
        raise NotImplementedError("Only able to load files of photomaps")
    
    return FlatFielder(photomaps, **kwargs)
