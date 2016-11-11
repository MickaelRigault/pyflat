#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" This module defines the basics objects used for the star flat analysis"""

import warnings
import numpy     as np

# - Astropy
from astropy import units

# - Astrobject
from astrobject.utils.decorators import _autogen_docstring_inheritance
from astrobject.baseobject import CatalogueHandler
from astrobject.collection import BaseCollection
from astrobject.collections.photospatial import get_sepobject, PhotoMap


def get_sepcollection( filenames, **kwargs):
    """ """
    list_photomaps = [get_sepobject(filename, **kwargs) for filename in filenames]    
    return PhotoMapCollection(list_photomaps, id_=np.range(len(filenames)))


class PhotoMapCollection( BaseCollection, CatalogueHandler ):
    """ Collection of PhotoMaps """

    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    def __init__(self, photomaps=None, empty=False, catalogue=None, **kwargs):
        """ """
        self.__build__()
        if empty:
            return

        # - Load the Catalogue if given
        if catalogue is not None:
            self.set_catalogue(catalogue, **kwargs)
            
        # - Load any given photomap
        if photomaps is not None:
            if hasattr(photomaps, "__iter__"):
                [self.add_photomap(photomap, **kwargs) for photomap in photomaps]
            else:
                self.add_photomap(photomaps, **kwargs)

    # =================== #
    #   Main Methods      #
    # =================== #
    def add_photomap(self, photomap, id_=None,
                     match_catalogue=True, **kwargs):
        """ 
        This method enables to load the given photomap in the
        self.photomaps container (dict).
        
        Parameters
        ----------
        new_image: [string or astrobject's Image (or children)]

        - option -
        id_: [any]
            key used to access the photomap from the _handler.
            id_ is set the photopoint's bandname if None

        **kwargs goes to catalogue matching.
        Return
        ------
        Void
        """
        # --------------
        # - Define ID
        if photomap is None:
            return

        if PhotoMap not in photomap.__class__.__mro__:
            raise TypeError("The givem photomap must be (or inherate of) an astrobject's PhotoMap")
        
            
        if self.has_sources():
            if id_ is None or id_ in self.list_id:
                id_= ""
                i = 1
                while id_+"%d"%i in self.list_id:
                    i+=1
                id_ = id_+"%d"%i
        else:
            id_ = "1"
            
        # -------------------
        # - Match the catalogue
        if self.has_catalogue() and match_catalogue:
            photomap.match_catalogue(force_it=True, **kwargs)
        # -------------------
        # - Record the map            
        self.photomaps[id_]  = photomap


    # -------------- #
    #  Catalogue     #
    # -------------- #
    def match_catalogue(self, deltadist=2*units.arcsec, **kwargs):
        """ Match the catalogue to all the PhotoMaps """
        if not self.has_catalogue():
            raise AttributeError("No catalogue to match")
        
        [self.photomaps[id_].match_catalogue(deltadist=2*units.arcsec,
                                              **kwargs)
         for id_ in self.list_id]


    def get_catindexes(self, inclusive=False):
        """ Get the index of the catalogue that has been mathced.
        if inclusive is False, only the catalogue entries matched in *all*
        individual photomaps will be returned. If False *any* catalogue
        entries that has been matched at least once will be returned

        Returns
        -------
        list (indexes of the catalogue)
        """
        catindexes = [self.photomaps[id_].catmatch["idx_catalogue"] for id_ in self.list_id]
        if inclusive:
            return np.unique(np.concatenate(catindexes))
        if len(catindexes) == 1:
            return catindexes[0]
        
        return list(frozenset(catindexes[0]).intersection(*catindexes[1:]))

    def catindex_to_index(self, catindex, dictformat=True, cleanindex=False):
        """ Get the index of the individual photomaps associated to the given catindex
        
        Parameters
        ----------
        catindex: [int of list of]
            catalogue entry that should be matched with the photomap indexes

        dictformat: [bool] -optional-
            The output will be returned in form of {id_: list_of_photomap_indexes}.
            If False this will be a list of list following list_id sorting.
            
        Returns
        -------
        dict following the ids
        """
        if dictformat:
            return {id_: self.photomaps[id_].catindex_to_index(catindex,cleanindex=cleanindex)
                    for id_ in self.list_id}
        else:
            return [self.photomaps[id_].catindex_to_index(catindex,cleanindex=cleanindex)
                    for id_ in self.list_id]

            
    # - Super Mother CatalogueHandler
    
    @_autogen_docstring_inheritance(CatalogueHandler.download_catalogue,"CatalogueHandler.download_catalogue")
    def download_catalogue(self, source="sdss",
                           set_it=True,force_it=False,
                           radec=None, radius_degree=None,
                           match_catalogue=True,match_angsep=2*units.arcsec,
                           **kwargs):
        #
        # default definition of radec and radius_degree
        #
        if radec is None or radius_degree is None:
            ra,dec = np.concatenate(self.get(["ra","dec"]), axis=0).T
            
        if radec is None:
            radec = np.mean([ra,dec],axis=1)
            radec = "%s %s"%(radec[0],radec[1])
        if radius_degree is None:
            radius_degree = np.max(np.std([ra,dec],axis=1)*5)
            
        out = super(PhotoMapCollection, self).download_catalogue(source=source,
                                                        set_it=set_it, force_it=force_it,
                                                        radec=radec, radius_degree=radius_degree,
                                                        **kwargs)
        if self.has_catalogue() and match_catalogue:
            self.match_catalogue(deltadist=match_angsep)

        return out
    @_autogen_docstring_inheritance(CatalogueHandler.set_catalogue,"CatalogueHandler.set_catalogue")
    def set_catalogue(self, catalogue, force_it=False, **kwargs):
        # 
        super(PhotoMapCollection, self).set_catalogue( catalogue, force_it=force_it, **kwargs)
        [self.photomaps[id_].set_catalogue(self.catalogue, force_it=force_it, **kwargs)
         for id_ in self.list_id]
        
    # =================== #
    #   Properties        #
    # =================== #
    @property
    def photomaps(self):
        """ """
        return self._handler
