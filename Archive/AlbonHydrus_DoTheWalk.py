# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:00:41 2016

@author: supersoil
"""

#==============================================================================
# Do the walk: Hydrus_DIY for all ALBON sites!
#==============================================================================

import numpy
from AlbonHydrus_SpecificFunk import *

#FC_L=['A','E','L','P','S','Z']
FC_L=['L','P','S','Z']
SFC_L=['REF']
HC_L=['H','C']
K_L=['KSK0','KS','K0']

for FC in FC_L:
    for SFC in SFC_L:
        for HC in HC_L:
            for K in K_L:
                print ('%s - %s - %s - %s'%(FC,SFC,HC,K))
                try: 
                    RunAlbonHydrus(FC,SFC,HC,K)
                except:
                    print('Unexpected error')

