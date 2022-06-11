# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 21:45:28 2020
@author: David Vasquez
Function: sellmeierDispForm
Optical Glass information took from:
    SCHOTT: Optical Glass - Datasheets - May 2019
    OHARA:  Glaskatalog_komplett - Version_19_Okt_2018
"""
try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))

import math
#Dictionary with Sellmeier Constants for Dispersion Formula
#/////////////NAME///////B1//////////B2//////////B3/////////C1////////////C2///////////C3///////////       
OPT_GLASS = { 
             'N-FK5'   :[0.84430934    ,0.344147824   ,0.910790213   ,0.00475111955 ,0.0149814849  , 97.8600293   ],
             'N-PSK3'  :[0.88727211    ,0.489592425   ,1.048652960   ,0.00469824067 ,0.0161818463  ,104.3749750   ],
             'N-BK7'   :[1.03961212    ,0.231792344   ,1.010469450   ,0.00600069867 ,0.0200179144  ,103.5606530   ],
             'N-K5'    :[1.08511833    ,0.199562005   ,0.930511663   ,0.00661099503 ,0.0241108660  ,111.9827770   ],
             'N-BAF10' :[1.58514950    ,0.143559385   ,1.085212690   ,0.00926681282 ,0.0424489805  ,105.6135730   ],
             'N-SF5'   :[1.52481889    ,0.187085527   ,1.427290150   ,0.01125475600 ,0.0588995392  ,129.1416750   ],
             'N-SF10'  :[1.62153902    ,0.256287842   ,1.644475520   ,0.01222414570 ,0.0595736775  ,147.4687930   ],
             'S-BAH11' :[1.57138860E+00,1.47869313E-01,1.28092846E+00,9.10807936E-03,4.02401684E-02,1.30399367E+02]
             }

def sellmeierDispForm (glass,wvln):
    '''
    glass:(list) List of the Sellmeier dispersion constants [B1.B2,B3,C1,C2,C3]
    wnln :(float)Wavelength in nm
    '''
    # Convert from nm to µm
    wvln2 = (wvln/1000.0)**2                 
    term1 = (glass[0]*wvln2)/(wvln2-glass[3])
    term2 = (glass[1]*wvln2)/(wvln2-glass[4])
    term3 = (glass[2]*wvln2)/(wvln2-glass[5])
    
    refIndex= math.sqrt(1 + term1 + term2 + term3)
    
    return refIndex

if __name__ == '__main__':
    coef= sellmeierDispForm (OPT_GLASS['N-SF5'],1060.0)